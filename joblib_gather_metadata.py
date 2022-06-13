#!/usr/bin/env python3
'''
Script to pull Uniprot Accession Identifier codes from blast sequence alignment results and then request and parse the metadata file associated with those Uniprot ID

USAGE: 
    python3 -gbff genbank_file -blast blast_output -cut 100 -spmn 9DELT -c 10

INPUT:
    gbff-file: string, path to the genbank gbff file
    blast_results: string, path to the blast output file
    cut: float, sets the sequence identity cutoff, which is applied to blast hits; below this cutoff and the balast hit is ignored
    spmn: string, used in uniprot searches for a specific species mnemonic
    c: number of threads to be made available

OUTPUT:
    Saves a pickle file of a list of dictionaries. each dictionary is filled with quantitative and metadata associated with a set of proteins
        [0]: uniprot_dictionary is a dictionary of dictionaries, keys based on WP codes, values are dictionaries filled with some of the metadata pulled from the uniprot metadata file
        [1]: proteinID_to_uniprotID is a dictionary of lists, keys based on WP codes, values are lists of lists, each element in the list is a list associated with a blast hit. 
        [2] and [3]: prot_to_loc and loc_to_prot are dictionaries giving the mappings between WP codes and locus tags (new and old). 
'''

import argparse
import datetime
import pickle
from urllib import request

from Bio import SeqIO, Blast, SearchIO

from joblib import Parallel, delayed

#######################################
### PARSING FUNCTIONS
#######################################

def uniprot_search(uniprotID,species_mnemonic=None,log_file_handle=None):
    '''code to access UniprotKB RESTful API metadata file for a uniprotID
    pull _some_ of the relevant metadata and store in a dictionary 
    '''
    metadata_url = 'https://www.uniprot.org/uniprot/%s.txt'
    dictionary_space = {'uniprotID': uniprotID,
                        'entry_name': '',
                        'status': '',
                        'sequence_length': '',
                        'original_date': '',
                        'sequence_date': '',
                        'modified_date': '',
                        'database_references': [],
                        'success': 0
                        }
    
    try:
        response = request.urlopen(metadata_url%(uniprotID))
        response_str = str(response.read())
        response_lines = response_str[2:-1].split('\\n')
        for line in response_lines:
            if line[:2] == 'ID':
                temp = line.split()
                if species_mnemonic and species_mnemonic not in temp[1]:
                    log_file_handle.write(f'This uniprotID ({uniprotID}) does not originate from the expected species ({species_mnemonic}). Killing task.\n')
                    break  # pops us out of the loop and returns an empty dictionary space with a 'success' value of 0...
                
                dictionary_space['success'] = 1
                dictionary_space['entry_name'] = temp[1]
                dictionary_space['status']     = temp[2][:-1]           # there's a stupid semicolon at the end of the status 
                dictionary_space['sequence_length'] = int(temp[3])
            elif line[:2] == 'DT':
                if 'integrated' in line:
                    dictionary_space['original_date'] = datetime.datetime.strptime(line.split()[1][:-1],'%d-%b-%Y').strftime('%Y-%m-%d')
                elif 'sequence version' in line:
                    dictionary_space['sequence_date'] = datetime.datetime.strptime(line.split()[1][:-1],'%d-%b-%Y').strftime('%Y-%m-%d')
                elif 'entry version' in line:
                    dictionary_space['modified_date'] = datetime.datetime.strptime(line.split()[1][:-1],'%d-%b-%Y').strftime('%Y-%m-%d')
            elif line[:2] == 'DR':
                dictionary_space['database_references'].append([line[5:]])
        return dictionary_space
    
    except Exception as e:
        e = str(e)
        if '400' in e:
            log_file_handle.write(uniprotID, e, 'Bad request. There is a problem with your input.')
        elif '404' in e:
            log_file_handle.write(uniprotID, e, "Not found. The resource you requested doesn't exist.")
        elif '410' in e:
            log_file_handle.write(uniprotID, e, 'Gone. The resource you requested was removed.')
        elif '500' in e:
            log_file_handle.write(uniprotID, e, 'Internal server error. Most likely a temporary problem, but if the problem persists please contact us.')
        elif '503' in e:
            log_file_handle.write(uniprotID, e, 'Service not available. The server is being updated, try again later.')
        else:
            log_file_handle.write(uniprotID, e, 'something else funky is going on...')
        return dictionary_space

def parse_blast(blast_file,seqid_test = False, seqid_cutoff = 90.0):
    """parse a blast output file with top hits for each query sequence
    gather quantitative metrics and metadata for each hit
    """
    proteinID_to_uniprotID = {}
    proteinID = 'hello world'
    with open(blast_file,'r') as blast_results:
        for line in blast_results:
            if "Query:" in line:
                temp = line.split()
                proteinID = temp[2]
                nRes = temp[3]
                proteinID_to_uniprotID[proteinID] = []
            elif proteinID in line:
                # check if seqid_test is being performed
                if seqid_test:
                    # perform the seqid test; greater than seqid cutoff
                    if float(line.split()[2]) >= seqid_cutoff:
                        proteinID_to_uniprotID[proteinID].append(line.split())
                    else:
                        continue
                else: 
                    proteinID_to_uniprotID[proteinID].append(line.split())
            else: 
                continue
    return proteinID_to_uniprotID

def protein_locus_dicts(genbank_file):
    """ This returns dicts for mapping between protein ID and locus tags from
        a Genbank file

    Note that we ignore the mapping from protein to locus tag for the old locus
    tags because then there would be multiple possible locus tags per protein.
    And, generally, we lookup locus tags for proteins, not the other way round.

    :param genbank_file: is the Genbank file from which we want to read
    :returns: two dicts, protein ID -> locus tag and locus tag -> protein ID
    """
    records = list(SeqIO.parse(genbank_file, "genbank"))

    locus_to_protein = {}
    protein_to_locus = {}

    for record in records:
        for feature in record.features:
            if 'protein_id' in feature.qualifiers:
                # testing version; the list will first be filled with the locus_tag feature and then the old_locus_tag feature
                protein_to_locus[feature.qualifiers['protein_id'][0]] = []
                
                # We have a protein ID, but there may be multiple incarnations
                # for a corresponding locus tag in this record.  We try
                # all the possible combos and create corresponding look-up
                # entries for every single one of them.  Oh, and we strip out
                # any locus tag underscores so that we have a consistent,
                # universal standard for matching locus tags with their proteins
                if 'locus_tag' in feature.qualifiers:
                    for locus_tag in feature.qualifiers['locus_tag']:
                        locus_to_protein[locus_tag.replace('_','')] = feature.qualifiers['protein_id'][0]
                        #protein_to_locus[feature.qualifiers['protein_id'][0]] = [locus_tag.replace('_','')]
                        protein_to_locus[feature.qualifiers['protein_id'][0]].append(locus_tag.replace('_',''))

                if 'old_locus_tag' in feature.qualifiers:
                    for old_locus_tag in feature.qualifiers['old_locus_tag']:
                        locus_to_protein[old_locus_tag.replace('_','')] = feature.qualifiers['protein_id'][0]
                        protein_to_locus[feature.qualifiers['protein_id'][0]].append(old_locus_tag.replace('_',''))
                        # This will *over-write* any previous mappings from
                        # the protein ID to locus tag, which we do not want.
                        # protein_to_locus[feature.qualifiers['protein_id'][0]] = \
                        #     old_locus_tag.replace('_','')

    return protein_to_locus, locus_to_protein

#######################################
### SUBMISSION FUNCTIONS
#######################################

def uniprot_search_pipeline(wp_code,blast_hits,species_mnemonic=''):
    """
    """

    if len(blast_hits) == 0:
        return [wp_code,{'success':0}] # no blast hits to be searched for... returning an empty dictionary
    else:
        for blast_hit in blast_hits:
            temp_dict = uniprot_search(blast_hit[1],species_mnemonic=species_mnemonic)
            if temp_dict['success']:
                return [wp_code,temp_dict]  # found the uniprot metadata for the respective blast hit
            else:
                print(wp_code, blast_hit[1], f'do not match up correctly based on expected species mnemonic ({species_mnemonic})')
                continue
        return [wp_code,{'success':0}] # never found an acceptable uniprot ID... returning an empty dictionary

#######################################
### MAIN
#######################################
if __name__ == '__main__':
    # read command line arguments.
    parser = argparse.ArgumentParser(description='Gather metadata information from a variety of files; store the results in a pkl file')
    parser.add_argument('--gbff-file', '-gbff', required=True, help='string, path to the genbank gbff file')
    parser.add_argument('--blast-results', '-blast', required=True, help='string, path to the blast output file')
    parser.add_argument('--seqID-cutoff', '-cut', default=0.0, type=float, help='float, sets the sequence identity cutoff, which is applied to blast hits; below this cutoff and the balast hit is ignored')
    parser.add_argument('--species-mnemonic', '-spmn', default=None, help="string, used in uniprot searches for a specific species mnemonic")
    parser.add_argument('--max-threads', '-c', required=True, type=int, help='number of threads to be made available')
    args = parser.parse_args()

    #parse genbank file
    prot_to_loc, loc_to_prot = protein_locus_dicts(args.gbff_file)
    #parse blast out file
    proteinID_to_uniprotID = parse_blast(args.blast_results,seqid_test=True,seqid_cutoff=args.seqID_cutoff)

    successes     = []
    no_blast_hits = []
    uniprot_fails = []

    log_file_handle = open('uniprot_search.log','w')
    results = Parallel(n_jobs=args.max_threads,prefer="threads")(
            delayed(uniprot_search_pipeline)(key,proteinID_to_uniprotID[key],species_mnemonic=args.species_mnemonic,log_file=log_file_handle) for key in proteinID_to_uniprotID.keys())

    # need to turn the list of lists of [key, dictionary] values into a dictionary of dictionaries
    no_blast_hits = []
    uniprot_fails = []
    successes     = []
    uniprot_dictionary = {}
    for result in results:
        wp_code = result[0]
        if result[1]['success']:
            uniprot_dictionary[wp_code] = result[1]
            successes.append(wp_code)
        elif len(proteinID_to_uniprotID[wp_code]) == 0:
            no_blast_hits.append(wp_code)
        else: 
            uniprot_fails.append(wp_code)

    log_file_handle.write( f'### {len(successes)} proteins successfully were parsed and IDed in the Uniprot.\n### {len(no_blast_hits)} proteins did not have any blast hits pass the sequence identity cutoff ({args.seqID_cutoff})\n### {len(uniprot_fails)} proteins did not have any successful uniprot hits whether this is due to query failures or if the uniprotIDs did not match the species mnemonic ({args.species_mnemonic})'))
    log_file_handle.close()

    today = datetime.date.today()
    with open('parsing_overview_%s.lst'%(today.strftime('%Y-%m-%d')),'w') as output:
        output.write('# successes\n')
        for key in successes:
            output.write(key+'\n')
        output.write('# no blast hits\n')
        for key in no_blast_hits:
            output.write(key+'\n')
        output.write('# uniprot fails\n')
        for key in uniprot_fails:
            output.write(key+'\n')

    with open('parsing_results.pkl','wb') as handle:
        pickle.dump([uniprot_dictionary,proteinID_to_uniprotID,prot_to_loc,loc_to_prot],handle)
        # saved as a list of dictionaries. 
        # uniprot_dictionary is a dictionary of dictionaries, keys based on WP codes, values are dictionaries filled with some of the metadata pulled from the uniprot metadata file
        # proteinID_to_uniprotID is a dictionary of lists, keys based on WP codes, values are lists of lists, each element in the list is a list associated with a blast hit. 
        # prot_to_loc and loc_to_prot are dictionaries giving the mappings between WP codes and locus tags (new and old). 

