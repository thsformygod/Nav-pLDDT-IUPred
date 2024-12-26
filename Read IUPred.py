# This script is basic to run IUPred.py on input sequence, and return IUPred values
# This need iupred3.py, please install IUPred at first, it could be found by https://iupred3.elte.hu/download_new




import tempfile
import os
import subprocess


def run_iupred3_long(input_sequence):
    # Create a temporary FASTA file
    with tempfile.NamedTemporaryFile(mode='w', delete=False) as temp_fasta:
        temp_fasta.write(">input_sequence\n")
        temp_fasta.write(input_sequence)

    # Run iupred3.py with the temporary FASTA file
    cmd = f"python D:\IUPred\iupred3.py {temp_fasta.name} long"
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()

    # Clean up the temporary FASTA file
    os.remove(temp_fasta.name)

    return stdout.decode(), stderr.decode()


def parse_iupred_output(output):
    list1, list2, list3 = [], [], []
    lines = output.strip().split('\n')

    # Parse data lines (skipping lines starting with "#")
    for line in lines:
        if not line.startswith("#"):
            parts = line.split("\t")
            if parts[0]:
                list1.append(int(parts[0]))
                list2.append(parts[1])
                list3.append(float(parts[2]))
            else:
                pass

    return list1, list2, list3

def process_input_file(sequence):
    stdout, stderr = run_iupred3_long(sequence)
    _, _, list3 = parse_iupred_output(stdout)
    return list3

# example protein sequence P00963
example_sequence = 'MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGCEKAVQVKVKALPDAQFEVVHSLAKWKRQTLGQHDFSAGEGLYTHMKALRPDEDRLSPLHSVYVDQWDWERVMGDGERQFSTLKSTVEAIWAGIKATEAAVSEEFGLAPFLPDQIHFVHSQELLSRYPDLDAKGRERAIAKDLGAVFLVGIGGKLSDGHRHDVRAPDYDDWSTPSELGHAGLNGDILVWNPVLEDAFELSSMGIRVDADTLKHQLALTGDEDRLELEWHQALLRGEMPQTIGGGIGQSRLTMLLLQLPHIGQVQCGVWPAAVRESVPSLL'

IUPred_score = process_input_file(example_sequence)
print(IUPred_score)