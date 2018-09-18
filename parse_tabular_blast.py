# This file contains a set of functions for parsing out some useful information
# from BLAST results files saved in BLAST's tabular output format ("-outfmt 6").

# Biopython is required for reading multifasta files and storing sequences.
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

# if all of your genome sequences are within one multifasta file
recs = [rec for rec in SeqIO.parse('all_genomes.fasta', 'fasta')]

# if they are within multiple multifasta files:
# fasta_files contains the paths to all of your multifasta files
fasta_files = [open('test_gg_97.fa', 'r')] 
recs = [rec for f in fasta_files for rec in SeqIO.parse(f, 'fasta')] 
# in this case, each multifasta could have multiple contigs so each contig may
# have a different id/name. This id/name would be it's identifier for BLAST.

# get a dictionary of your strain/contig names and the corresponding SeqRecord 
# so that you can retrieve the SeqRecord using the name of that SeqRecord
strain_name_seq_rec = {rec.name:rec for rec in recs}


# The lengths can be stored in a dictionary where the query id/name can be
# used to retrieve the length of that query.
# query_files contains all of the query multifasta file paths
query_files = [open('16s_rrna_db.faa', 'r')]
query_lengths = {rec.name:len(rec.seq) for q_file in query_files for rec in SeqIO.parse(q_file, 'fasta')}

def adjust_subject_indices(q_s_i, q_e_i, s_s_i, s_e_i, q_len, revcomp):
    """
    Adjust subject indices based on query indices and query length.

    Args:
        q_s_i (int): Query start index
        q_e_i (int): Query end index
        s_s_i (int): Subject start index
        s_e_i (int): Subject end index
        q_len (int): Query length
        revcomp (bool): Reverse complemented?

    Returns:
        A tuple containing the adjusted subject start and end indices.
    """
    add_to_start = q_s_i - 1
    add_to_end = q_len - q_e_i
    if revcomp:
        s_s_i -= add_to_end
        s_e_i += add_to_start
    else:
        s_s_i -= add_to_start
        s_e_i += add_to_end
    s_s_i -= 1
    s_e_i -= 1
    return (s_s_i, s_e_i)

def extract_extended_seq(sseq, q_s_i, q_e_i, s_s_i, s_e_i, q_len):
    """
    Extract the extended subject sequence based on query start and end indices.

    Args:
        sseq (Bio.Sequence,string): Subject sequence
        q_s_i (int): Query start index
        q_e_i (int): Query end index
        s_s_i (int): Subject start index
        s_e_i (int): Subject end index
        q_len (int): Query length

    Returns:
        The extended subject sequence that matches up with the query 
    """
    revcomp = s_s_i > s_e_i
    s_s_i, s_e_i = adjust_subject_indices(
        q_s_i, 
        q_e_i, 
        s_s_i, 
        s_e_i, 
        q_len, 
        revcomp)
    
    if s_s_i < 0:
        if len(sseq) < s_e_i:
            return sseq
        else:
            return sseq[0: s_e_i + 1]
    else:
        if len(sseq) - s_s_i < s_e_i - s_s_i + 1:
            return sseq[s_s_i:len(sseq) - s_s_i]
        else:
            return sseq[s_s_i:s_e_i + 1]
    

def parse_tabular_blast_results(blast_results_file):
    """
    Parse out the BLAST results from a tabular format BLAST results file into a 
    2D dictionary containing query id/name -> subject id/name -> BLAST result.

    Args:
        blast_results_file (string): Path to tabular BLAST results file.

    Returns:
        2D dictionary containing query id/name -> subject id/name -> BLAST result
    """
    # results are stored in a 2D dictionary - query -> subject -> BLAST result
    results = {}
    with open(blast_results_file,'r') as f:
        for line in f:
            line = line.rstrip()
            sp = line.split('\t')
            # get query id
            q = sp[0]
            q_len = query_lengths[q]
            # get subject id
            s = sp[1]
            # get the subject SeqRecord object
            sseq = strain_name_seq_rec[s]
            # get % identity and convert to decimal (e.g. 96.7 -> 0.967)
            pid = float(sp[2])/100.0
            # get alignment length
            aln_len = int(sp[3])
            # get query start index
            q_s_i = int(sp[6])
            # get query end index
            q_e_i = int(sp[7])
            # get subject start index
            s_s_i = int(sp[8])
            # get subject end index
            s_e_i = int(sp[9])
            # try to get the extended subject sequence based on the full length 
            # of the query. If not able to grab the full length, grab as much
            # as possible. 
            extracted_seq = extract_extended_seq(sseq, 
                q_s_i, 
                q_e_i, 
                s_s_i, 
                s_e_i,
                q_len)
            # number of identities (count of positive sequence matches)
            identities = int(aln_len*pid)
            # the results for a particular query could be whatever you need
            # in this example, a dictionary is made where different results have
            # different identifiers for instance
            result_dict = {
            'aln_len':aln_len, 
            'q_len':q_len,
            'pid':pid, 
            'extracted_seq': extracted_seq
            }
            # could add more items to this BLAST result dictionary
            if q not in results:
                results[q] = {s:[result_dict]}
            else:
                if s not in results[q]:
                    results[q][s] = [result_dict]
                else:
                    results[q][s].append(result_dict)
    return results

# Parsing the results file from BLASTN.
parse_tabular_blast_results(open('test_out.txt', 'r'))
