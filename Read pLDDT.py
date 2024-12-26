# This script is basicly to read AF files and return sequence and pLDDT scores

import os

def read_AF_cif(input_file):
    if os.path.exists(input_file):
        filename = os.path.basename(input_file)
        pLDDT = []
        AF_sequence = []
        amino_acids_dict = {'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys', 'Q': 'Gln', 'E': 'Glu',
                            'G': 'Gly', 'H': 'His', 'I': 'Ile', 'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe',
                            'P': 'Pro', 'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val', 'X': 'Xaa'}

        amino_acids_dict = {key: value.upper() for key, value in amino_acids_dict.items()}

        amino_acids_dict = {value: key for key, value in amino_acids_dict.items()}

        # Read the CIF file
        with (open(input_file, 'r') as cif_content):
            lines = cif_content.readlines()
            inside_data_loop = False
            # Process each line in the CIF file
            for line in lines:
                # Check if the line contains "_ma_qa_metric_local.ordinal_id"
                if "_ma_qa_metric_local.ordinal_id" in line:
                    inside_data_loop = True
                    continue
                # Check if we are inside the relevant data loop
                if inside_data_loop:
                    # Check if the line is empty or starts with "#", indicating the end of the data loop
                    if not line.strip() or line.startswith("#"):
                        break
                    # Check if the line contains label_comp_id and metric_value information
                    columns = line.split()
                    #print(columns)
                    label_comp_id = columns[3]
                    metric_value = round(float(columns[4]),1)
                    pLDDT.append(metric_value)

                    AA = columns[1].upper()
                    #print ("current AA is " + AA)
                    AA = amino_acids_dict.get(AA, "X")

                    #print("dicted AA is " + AA)
                    AF_sequence.append(AA)


            AF_sequence = ''.join(str(num) for num in AF_sequence)

            return pLDDT,AF_sequence
    else:
        error_message = "doesnt_exit"
        return error_message,error_message

input_example = "D:\\Examples\\02_read_extract\\example\\inputs\\AF-P43308-F1-model_v4.cif"

pLDDT,sequence = read_AF_cif(input_example)
print(pLDDT)
print(sequence)
