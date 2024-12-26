# original script name: read information from cif03_02.py

# This script is used to read PDB .cif files, and extract related information, such as R-factors, Resolution, id of entries, reference ids and so on.

#version: biopython 1.79


# How to use:
# input_file01 = 'path\PDB_ids.txt'
# This PDB_ids.txt should contain one column of PDBs. as example_PDB_ids.txt


# input_path01 = 'path_PDBs\'
# path_PDBs is where put unziped PDBs, to let read_cit function go and read PDBs

# output_file = "path\output.txt"


from Bio.PDB.MMCIF2Dict import MMCIF2Dict
import os
def read_cif(cif_filename):
    with open(cif_filename, 'r', encoding='utf-8') as handle:
        mmcif_dict = MMCIF2Dict(handle)
        entry_id_value = '|'.join(mmcif_dict["_entry.id"])if "_entry.id" in mmcif_dict else "None"
        entry_id_value = "xxxx"+ str(entry_id_value)+"xxxx"
        method = '|'.join(mmcif_dict["_exptl.method"]) if "_exptl.method" in mmcif_dict else "None"
        ls_d_res_low = '|'.join(mmcif_dict["_refine.ls_d_res_low"]) if "_refine.ls_d_res_low" in mmcif_dict else "None"
        ls_d_res_high = '|'.join(mmcif_dict["_refine.ls_d_res_high"]) if "_refine.ls_d_res_high" in mmcif_dict else "None"
        ls_R_factor_obs = '|'.join(mmcif_dict["_refine.ls_R_factor_obs"]) if "_refine.ls_R_factor_obs" in mmcif_dict else "None"
        ls_R_factor_R_free = '|'.join(mmcif_dict["_refine.ls_R_factor_R_free"]) if "_refine.ls_R_factor_R_free" in mmcif_dict else "None"
        ls_R_factor_R_work = '|'.join(mmcif_dict["_refine.ls_R_factor_R_work"]) if "_refine.ls_R_factor_R_work" in mmcif_dict else "None"
        em_3d_reconstruction_resolution = '|'.join(mmcif_dict["_em_3d_reconstruction.resolution"]) if "_em_3d_reconstruction.resolution" in mmcif_dict else "None"
        em_experiment_reconstruction_method = '|'.join(mmcif_dict["_em_experiment.reconstruction_method"]) if "_em_experiment.reconstruction_method" in mmcif_dict else "None"

        recvd_initial_deposition_date = '|'.join(mmcif_dict["_pdbx_database_status.recvd_initial_deposition_date"])if "_pdbx_database_status.recvd_initial_deposition_date" in mmcif_dict else "None"
        release_date = '|'.join(mmcif_dict["_pdbx_audit_revision_history.revision_date"])if "_pdbx_audit_revision_history.revision_date" in mmcif_dict else "None"
        #title = '|'.join(mmcif_dict["_citation.title"])if "_citation.title" in mmcif_dict else "None"
        year = '|'.join(mmcif_dict["_citation.year"])if "_citation.year" in mmcif_dict else "None"
        journal_abbrev = '|'.join(mmcif_dict["_citation.journal_abbrev"])if "_citation.journal_abbrev" in mmcif_dict else "None"
        struct_ref_entity_id = '|'.join(mmcif_dict["_struct_ref.entity_id"])if "_struct_ref.entity_id" in mmcif_dict else "None"
        struct_ref_pdbx_db_accession = '|'.join(mmcif_dict["_struct_ref.pdbx_db_accession"]) if "_struct_ref.pdbx_db_accession" in mmcif_dict else "None"

        struct_ref_pdbx_db_accession = "x["+ str(struct_ref_pdbx_db_accession) + "]"
        align_id = '|'.join(mmcif_dict["_struct_ref_seq.align_id"]) if "_struct_ref_seq.align_id" in mmcif_dict else "None"
        pdbx_PDB_id_code = '|'.join(mmcif_dict["_struct_ref_seq.pdbx_PDB_id_code"]) if "_struct_ref_seq.pdbx_PDB_id_code" in mmcif_dict else "None"
        seq_align_beg = '|'.join(mmcif_dict["_struct_ref_seq.seq_align_beg"])if "_struct_ref_seq.seq_align_beg" in mmcif_dict else "None"
        seq_align_end = '|'.join(mmcif_dict["_struct_ref_seq.seq_align_end"])if "_struct_ref_seq.seq_align_end" in mmcif_dict else "None"
        pdbx_db_accession = '|'.join(mmcif_dict["_struct_ref_seq.pdbx_db_accession"])if "_struct_ref_seq.pdbx_db_accession" in mmcif_dict else "None"
        db_align_beg = '|'.join(mmcif_dict["_struct_ref_seq.db_align_beg"])if "_struct_ref_seq.db_align_beg" in mmcif_dict else "None"
        db_align_end = '|'.join(mmcif_dict["_struct_ref_seq.db_align_end"])if "_struct_ref_seq.db_align_end" in mmcif_dict else "None"


        ref_id = '|'.join(mmcif_dict["_struct_ref_seq.ref_id"]) if "_struct_ref_seq.ref_id" in mmcif_dict else "None"
        pdbx_strand_id = '|'.join(mmcif_dict["_struct_ref_seq.pdbx_strand_id"]) if "_struct_ref_seq.pdbx_strand_id" in mmcif_dict else "None"

        ################################################################
        entity_poly_entity_id = '|'.join(mmcif_dict["_entity_poly.entity_id"]) if "_entity_poly.entity_id" in mmcif_dict else "None"
        entity_poly_pdbx_strand_id = '|'.join(mmcif_dict["_entity_poly.pdbx_strand_id"]) if "_entity_poly.pdbx_strand_id" in mmcif_dict else "None"





        struct_ref_seq_pdbx_strand_id= '|'.join(mmcif_dict["_struct_ref_seq.pdbx_strand_id"]) if "_struct_ref_seq.pdbx_strand_id" in mmcif_dict else "None"
        struct_ref_seq_pdbx_PDB_id_code= '|'.join(mmcif_dict["_struct_ref_seq.pdbx_PDB_id_code"]) if "_struct_ref_seq.pdbx_PDB_id_code" in mmcif_dict else "None"


        seq_align_beg = '|'.join(mmcif_dict["_struct_ref_seq.seq_align_beg"])if "_struct_ref_seq.seq_align_beg" in mmcif_dict else "None"
        seq_align_end = '|'.join(mmcif_dict["_struct_ref_seq.seq_align_end"])if "_struct_ref_seq.seq_align_end" in mmcif_dict else "None"

        struct_ref_seq_pdbx_db_accession= '|'.join(mmcif_dict["_struct_ref_seq.pdbx_db_accession"])if "_struct_ref_seq.pdbx_db_accession" in mmcif_dict else "None"
        db_align_beg = '|'.join(mmcif_dict["_struct_ref_seq.db_align_beg"])if "_struct_ref_seq.db_align_beg" in mmcif_dict else "None"
        db_align_end = '|'.join(mmcif_dict["_struct_ref_seq.db_align_end"])if "_struct_ref_seq.db_align_end" in mmcif_dict else "None"



        return (entry_id_value,method,ls_d_res_low,ls_d_res_high,ls_R_factor_obs,ls_R_factor_R_free,ls_R_factor_R_work,em_3d_reconstruction_resolution,em_experiment_reconstruction_method,recvd_initial_deposition_date,release_date,year,journal_abbrev,struct_ref_entity_id,struct_ref_pdbx_db_accession,entity_poly_entity_id,entity_poly_pdbx_strand_id,struct_ref_seq_pdbx_strand_id,struct_ref_seq_pdbx_PDB_id_code,seq_align_beg,seq_align_end,struct_ref_seq_pdbx_db_accession,db_align_beg,db_align_end)


def read_cif_again(input_file01, input_path01, output_file):
    with (open(input_file01, 'r') as input_file01):
        for line in input_file01:
            line = line.strip("\n").split("\t")
            PDB = line[0]
            PDB_cif = input_path01 + PDB + ".cif"

            print(PDB)
            entry_id_value,method,ls_d_res_low,ls_d_res_high,ls_R_factor_obs,ls_R_factor_R_free,ls_R_factor_R_work,em_3d_reconstruction_resolution,em_experiment_reconstruction_method,recvd_initial_deposition_date,release_date,year,journal_abbrev,struct_ref_entity_id,struct_ref_pdbx_db_accession,entity_poly_entity_id,entity_poly_pdbx_strand_id,struct_ref_seq_pdbx_strand_id,struct_ref_seq_pdbx_PDB_id_code,seq_align_beg,seq_align_end,struct_ref_seq_pdbx_db_accession,db_align_beg,db_align_end = read_cif(
                PDB_cif)
            # Write the results to the output file
            with open(output_file, 'a') as output:
                output.write(f"{PDB}\t{entry_id_value}\t{method}\t{ls_d_res_low}\t{ls_d_res_high}\t{ls_R_factor_obs}\t{ls_R_factor_R_free}\t{ls_R_factor_R_work}\t{em_3d_reconstruction_resolution}\t{em_experiment_reconstruction_method}\t{recvd_initial_deposition_date}\t{release_date}\t{year}\t{journal_abbrev}\t{struct_ref_entity_id}\t{struct_ref_pdbx_db_accession}\t{entity_poly_entity_id}\t{entity_poly_pdbx_strand_id}\t{struct_ref_seq_pdbx_strand_id}\t{struct_ref_seq_pdbx_PDB_id_code}\t{seq_align_beg}\t{seq_align_end}\t{struct_ref_seq_pdbx_db_accession}\t{db_align_beg}\t{db_align_end}\n")

input_file01 = "D:\\Examples\\02_read_extract\\example\\inputs\\example_PDB_ids.txt"
input_path01 = "D:\\Examples\\02_read_extract\\example\\inputs\\"
output_file = "D:\\Examples\\02_read_extract\\example\\outputs\\example_PDB_ids_output.txt"

read_cif_again(input_file01, input_path01, output_file)