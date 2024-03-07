import rdflib as rl
import pandas as pd
from rdflib.namespace import RDF

## load ontology
g = rl.Graph()
g.parse("ImPO_new.owl",format="xml")
impo = "https://raw.githubusercontent.com/liseda-lab/ImPO/main/ImPO.owl#"
placeholder = rl.Namespace("http://placeholder.com/")

#######################################################
    ##################################################
    # transcriptomics DATA FILES INTEGRATION         #
    ##################################################
#######################################################


## Christophe data files (transcriptomics)
subjid = pd.read_csv("Christophe/SUBJID.csv", sep=",", index_col=None) # DONE
cell_deconv = pd.read_csv("Christophe/cell_deconvolution.csv", sep=",", index_col=None) # DONE
gene_mutations = pd.read_csv("Christophe/gene_mutations.csv", sep=",", index_col=None) # 
struc_var = pd.read_csv("Christophe/structural_variant.csv", sep=",", index_col=None) # DONE
pronostic = pd.read_csv("Christophe/pronostic_score.csv", sep=",", index_col=None) # DONE
tumor_feat = pd.read_csv("Christophe/tumor_feature.csv", sep=",", index_col=None) # 
bio_feat = pd.read_csv("Christophe/biological_feature.csv", sep=",", index_col=None) # DONE
bio_pathway = pd.read_csv("Christophe/biological_pathway.csv", sep=",", index_col=None) # DONE
gene_exp = pd.read_csv("Christophe/gene_expression.csv", sep=",", index_col=None) # DONE
drug_cat = pd.read_csv("Christophe/drug_category_and_molecule.csv", sep=",", index_col=None) # DONE
meta = pd.read_csv("Christophe/META_data.csv", sep=",", index_col=None) # DONE
hla_subtype = pd.read_csv("Christophe/HLA_subtype.csv", sep=",", index_col=None) # DONE

########################################################################################################################################
########################################################################################################################################
#
#
#
#
####################################### PATIENT

## SUBJID file, per patient
# for k in range(len(subjid["SUBJID"])):
#     patient_uri = rl.URIRef(placeholder + subjid["SUBJID"][k]) ## patient individuals uris
#     g.add((patient_uri, RDF.type, rl.URIRef(impo + "patient"))) # adds as individual of class patient
#     g.add((patient_uri, rl.URIRef(impo + "has_id"), rl.Literal(subjid["SUBJID"][k]))) # adds SUBJID
#     g.add((patient_uri, rl.URIRef(impo + "has_mskcc"), rl.Literal(subjid["MSKCC"][k]))) # adds MSKCC
#     g.add((patient_uri, rl.URIRef(impo + "has_imdc"), rl.Literal(subjid["IMDC"][k]))) # adds IMDC
#     g.add((patient_uri, rl.URIRef(impo + "has_sex"), rl.Literal(subjid["Sex"][k]))) # adds Sex
#     g.add((patient_uri, rl.URIRef(impo + "has_age"), rl.Literal(subjid["Age"][k]))) # adds Age
#     g.add((patient_uri, rl.URIRef(impo + "has_cohort"), rl.Literal(subjid["Cohort"][k]))) # adds Cohort
#     g.add((patient_uri, rl.URIRef(impo + "received_prior_therapy"), rl.Literal(subjid["Received_Prior_Therapy"][k])))  # adds Received_Prior_Therapy
#     g.add((patient_uri, rl.URIRef(impo + "has_number_of_prior_therapies"), rl.Literal(subjid["Number_of_Prior_Therapies"][k]))) # adds Number_of_Prior_Therapies
#     g.add((patient_uri, rl.URIRef(impo + "has_tmb_counts"), rl.Literal(subjid["TMB_Counts"][k]))) # adds TMB counts

# g.serialize(destination="transcriptomics_impo_subject.owl", format='xml') # generates owl file
# del g
# g = rl.Graph()
# g.parse("transcriptomics_impo_subject.owl")
# print("  ------  over 1")


# ##########################################################
# #
# #
# #
# #
# ######################################### SAMPLE
    
# # ## cell_deconvolution file,  per sample
# for k in range(len(cell_deconv["Sample_ID"])): # add sample individual uris
#     sample_uri = rl.URIRef(placeholder + cell_deconv["Sample_ID"][k].replace(" ","_")) ## sample individuals uris
#     g.add((sample_uri, RDF.type, rl.URIRef(impo + "sample"))) # adds as individual of class sample
#     g.add((sample_uri, rl.URIRef(impo + "has_id"), rl.Literal(cell_deconv["Sample_ID"][k]))) # adds sample_ID
#     g.add((sample_uri, rl.URIRef(impo + "has_pvalue"), rl.Literal(cell_deconv["P-value"][k]))) # adds p-value
#     g.add((sample_uri, rl.URIRef(impo + "has_correlation"), rl.Literal(cell_deconv["Correlation"][k]))) # adds sample_ID
#     g.add((sample_uri, rl.URIRef(impo + "has_rmse"), rl.Literal(cell_deconv["RMSE"][k]))) # adds p-value

# for k in range(len(meta["Sample_ID"])): # for samples in meta_data file:
#     sample_uri = rl.URIRef(placeholder + meta["Sample_ID"][k].replace(" ","_")) ## sample individuals uris
#     g.add((sample_uri, rl.URIRef(impo + "specimen_of"), rl.URIRef(placeholder + meta["SUBJID"][k]))) # adds relation "specimen_of" patient
#     sample_type_uri = rl.URIRef(placeholder + meta["Sample_Type"][k].replace(" ","_")) ## sample type individual uris
#     g.add((sample_uri, rl.URIRef(impo + "has_type"), sample_type_uri)) # links sample to sample type class

# g.serialize(destination="transcriptomics_impo_sample.owl", format='xml') # generates owl file
# # del g
# # g = rl.Graph()
# # g.parse("transcriptomics_impo_sample.owl")
# print("  ------  over 2")


import math
# subj_to_samp = {} ## THROUGH SAMPLE LINK WITH SUBJID THEN TO HLA_SUBTYPE TO CREATE HLA AND ADD ITS PROPERTIES IN SUBJID.CSV FILE (HLA is feature of some SAMPLE)
# for k in range(len(gene_mutations["Sample_ID"])): # get all samples
#     if type(gene_mutations["Sample_ID"][k]) != float or math.isnan(gene_mutations["Sample_ID"][k]) == False:
#         if subj_to_samp.get(gene_mutations["SUBJID"][k]) == None: # make dictionary with all samples per patient / SUBJID
#             subj_to_samp[gene_mutations["SUBJID"][k]] = []
#         if gene_mutations["Sample_ID"][k] not in subj_to_samp[gene_mutations["SUBJID"][k]]:
#             subj_to_samp[gene_mutations["SUBJID"][k]].append(gene_mutations["Sample_ID"][k])
# for k in range(len(meta["SUBJID"])):
#     if subj_to_samp.get(meta["SUBJID"][k]) == None:
#         subj_to_samp[meta["SUBJID"][k]] = []
#     if meta["Sample_ID"][k] not in subj_to_samp[meta["SUBJID"][k]]:
#         subj_to_samp[meta["SUBJID"][k]].append(meta["Sample_ID"][k])
# for k in range(len(hla_subtype["SUBJID"])):
#     if subj_to_samp.get(hla_subtype["SUBJID"][k]) != None: # for samples with HLA info
#         for j in subj_to_samp[hla_subtype["SUBJID"][k]]:
#             sample_uri = rl.URIRef(placeholder + j.replace(" ","_")) ## sample individuals uris
#             hla_uri = rl.URIRef(placeholder + hla_subtype["HLA_subtype"][k].split("*")[0]) ## HLA individuals uris       ### REVISIT THIS
#             g.add((hla_uri, RDF.type, rl.URIRef(impo + "HLA_" + hla_subtype["HLA_subtype"][k].split("*")[0]))) # adds as individual of class gene mutation
#             g.add((sample_uri, rl.URIRef(impo + "has_feature"), rl.Literal(hla_uri))) # adds HLA class
#             g.add((hla_uri, rl.URIRef(impo + "has_allele"), rl.Literal(hla_subtype["HLA_subtype"][k].split("*")[1]))) # adds allele property to hla

# g.serialize(destination="transcriptomics_impo_hla.owl", format='xml') # generates owl file
# # del g
# # g = rl.Graph()
# # g.parse("transcriptomics_impo_pre-4_step_5000.owl")
# print("  ------  over 3")


##########################################################
#
#
#
#
##################################################### GENE  AND  RNA-SEQ - GENE Quantification 
# ## gene_expression file, per gene 
# for k in range(len(gene_exp["gene_id"])): # add gene individual uris
#     gene_uri = rl.URIRef(placeholder + gene_exp["gene_id"][k]) ## gene individuals uris
#     g.add((gene_uri, RDF.type, rl.URIRef(impo + "gene"))) # adds as individual of class gene
#     g.add((gene_uri, rl.URIRef(impo + "has_id"), rl.Literal(gene_exp["gene_id"][k]))) # adds gene_id
#     g.add((gene_uri, rl.URIRef(impo + "has_name"), rl.Literal(gene_exp["hgnc_symbol"][k]))) # adds gene_id
#     for j in gene_exp.columns:
#         rna_seq_uri = rl.URIRef(placeholder + "assay_" + gene_exp["gene_id"][k] + "_" + j)
#         g.add((rna_seq_uri, RDF.type, rl.URIRef(impo + "RNA_sequencing"))) # adds as individual of class RNA sequencing
#         rna_seq_gene_quantification = rl.URIRef(placeholder + "rnaseq_gene_id_" + gene_exp["gene_id"][k] + "_" + j)
#         g.add((rna_seq_gene_quantification, RDF.type, rl.URIRef(impo + "RNA_sequencing-gene_quantification")))
#         g.add((rna_seq_gene_quantification, rl.URIRef(impo + "has_source"), rna_seq_uri))
#         g.add((rna_seq_gene_quantification, rl.URIRef(impo + "has_target"), gene_uri))
#         g.add((rna_seq_gene_quantification, rl.URIRef(impo + "has_count"), rl.Literal(gene_exp[j][k])))

#     if k % 10000 == 0:
#         g.serialize(destination="transcriptomics_gene_exp_{}.owl".format(k), format='xml')
#         g = rl.Graph()
#         g.parse("ImPO_new.owl",format="xml")
#         print("  ------  pre - 4 --> on step {}".format(k))

#     if k == (len(gene_exp["gene_id"]) - 1):
#         g.serialize(destination="transcriptomics_gene_exp_{}.owl".format(k), format='xml')


# # # g.parse("transcriptomics_impo_4.owl")
# print("  ------  over 4")



##########################################################
#
#
#
#
############################################ GENOMIC MUTATION

# ## STRUCTURAL_VARIANT file,  per genomic_mutation /CNV
# for k in range(len(struc_var["CNV_ID"])): # add genomic mutation individual uris
#     genomic_mutation_uri = rl.URIRef(placeholder + "genomic_mutation_" + struc_var["CNV_ID"][k]) ## genomic mutation individuals uris
#     g.add((genomic_mutation_uri, RDF.type, rl.URIRef(impo + "genomic_mutation_"))) # adds as individual of class genomic mutation
#     g.add((genomic_mutation_uri, rl.URIRef(impo + "has_id"), rl.Literal(struc_var["CNV_ID"][k]))) # adds CNV_ID
#     g.add((genomic_mutation_uri, rl.URIRef(impo + "has_amplification_peak"), rl.Literal(''.join([str(struc_var[j][k])+"," for j in struc_var.columns if "Amplification" in j and "CN" not in j]))))  # adds string of all amplification peak values
#     g.add((genomic_mutation_uri, rl.URIRef(impo + "has_deletion_peak"), rl.Literal(''.join([str(struc_var[j][k])+"," for j in struc_var.columns if "Deletion" in j and "CN" not in j]))))  # adds string of all deletion peak values
#     g.add((genomic_mutation_uri, rl.URIRef(impo + "has_amplification_peak-CN_values"), rl.Literal(''.join([str(struc_var[j][k])+"," for j in struc_var.columns if "Amplification" in j and "CN" in j]))))  # adds string of all amplification peak CN values
#     g.add((genomic_mutation_uri, rl.URIRef(impo + "has_deletion_peak-CN_values"), rl.Literal(''.join([str(struc_var[j][k])+"," for j in struc_var.columns if "Deletion" in j and "CN" in j]))))  # adds string of all deletion peak CN values

# g.serialize(destination="transcriptomics_genomic_mutation.owl", format='xml') # generates owl file
# # del g
# # g = rl.Graph()
# # g.parse("transcriptomics_impo_5.owl")
# print("  ------  over 5")


# #### HOW TO ADD PROPERTIES FROM GENE_MUTATIONS.CSV FILE TO GENE_MUTATION INDIVIDUALS:
#     ### CNV-ID IS A SAMPLE ID, GET SUBJID, IF MULTIPLE SAMPLE IDS FOR SUBJID, ADD CNV_ID PER SAMPLE, IF NO SAMPLE ID THEN ADD PER SUBJID TO ADD PROPERTIES TO CNV ID
# cnv_ids = {}
# for k in range(len(struc_var["CNV_ID"])):
#     cnv_ids[struc_var["CNV_ID"][k]] = []
# samp_to_subj = {}
# for k in range(len(gene_mutations["SUBJID"])): # get all samples
#     if type(gene_mutations["SUBJID"][k]) != float or math.isnan(gene_mutations["SUBJID"][k]) == False:
#         if samp_to_subj.get(gene_mutations["Sample_ID"][k]) == None: # make dictionary with all samples per patient / SUBJID
#             samp_to_subj[gene_mutations["Sample_ID"][k]] = []
#         if gene_mutations["SUBJID"][k] not in samp_to_subj[gene_mutations["Sample_ID"][k]]:
#             samp_to_subj[gene_mutations["Sample_ID"][k]].append(gene_mutations["SUBJID"][k])
# for k in range(len(meta["SUBJID"])):
#     if samp_to_subj.get(meta["Sample_ID"][k]) == None:
#         samp_to_subj[meta["Sample_ID"][k]] = []
#     if meta["SUBJID"][k] not in samp_to_subj[meta["Sample_ID"][k]]:
#         samp_to_subj[meta["Sample_ID"][k]].append(meta["SUBJID"][k])
# subj_dict = {}
# for k in range(len(meta["SUBJID"])):                                 # CNV_ID = SAMPLE ID from META_data.csv
#     if cnv_ids.get(meta["Sample_ID"][k]) != None:                    # if cnv_id is in samples,
#         if samp_to_subj.get(meta["Sample_ID"][k]) != None:           # if SUBJID linked to sample/cnv_id is in gene_mutations.csv file
#             for l in samp_to_subj[meta["Sample_ID"][k]]:
#                 if subj_dict.get(l) == None: subj_dict[l] = []     
#                 subj_dict[l].append(meta["Sample_ID"][k])            # save SUBJIDs linked to the cnv in a dictionary
# del samp_to_subj
# # del subj_to_samp
name_to_id = {}
for k in range(len(gene_exp["gene_id"])):
    if type(gene_exp["hgnc_symbol"][k]) != str:
        if math.isnan(gene_exp["hgnc_symbol"][k]) != True:
            if name_to_id.get(gene_exp["hgnc_symbol"][k]) == None: name_to_id[gene_exp["hgnc_symbol"][k]] = []
            name_to_id[gene_exp["hgnc_symbol"][k]].append(gene_exp["gene_id"][k])

# #### INSTANCIATE PROTEIN and GENE
# for k in range(len(gene_mutations["gene_name"])):
#     try:
#         gene_uri = rl.URIRef(placeholder + name_to_id[gene_mutations["gene_name"][k]][0])
#         g.add((gene_uri, RDF.type, rl.URIRef(impo + "gene")))
#         g.add((gene_uri, rl.URIRef(impo + "has_id"), rl.Literal(name_to_id[gene_mutations["gene_name"][k]][0])))
#         g.add((gene_uri, rl.URIRef(impo + "has_name"), rl.Literal(gene_mutations["gene_name"][k])))
#         protein_uri = rl.URIRef(placeholder +  "gene_product_" + name_to_id[gene_mutations["gene_name"][k]][0])
#         g.add((protein_uri, RDF.type, rl.URIRef(impo + "protein")))
#         # g.add((protein_uri, rl.URIRef(impo + "has_id"), rl.Literal(gene_mutations["protein id"][k])))
#     except: pass

# g.serialize(destination="transcriptomics_impo_gene_prot.owl", format='xml') # generates owl file
# # del g
# # g = rl.Graph()
# # g.parse("transcriptomics_impo_6.owl")
# print("  ------  over 6")


# for k in range(len(gene_mutations["SUBJID"])): ## MAKE AND INSTANTIATE CLASS GENOMIC MUTATION
#     if subj_dict.get(gene_mutations["SUBJID"][k]) != None: # get correct CNV_ID
#         for m in subj_dict[gene_mutations["SUBJID"][k]]:
#             cnv_id_uri = rl.URIRef(placeholder + "genomic_mutation_" + m)
#             g.add((cnv_id_uri, rl.URIRef(impo + "has_reference_allele"), rl.Literal(gene_mutations["Reference_Allele"][k]))) # adds reference allele to genomic mutation
#             g.add((cnv_id_uri, rl.URIRef(impo + "has_chromosome_number"), rl.Literal(gene_mutations["Chromosome"][k]))) # adds chromosome number
#             g.add((cnv_id_uri, rl.URIRef(impo + "has_codon_change"), rl.Literal(gene_mutations["Codon_Change"][k])))
#             g.add((cnv_id_uri, rl.URIRef(impo + "has_cDNA_change"), rl.Literal(gene_mutations["cDNA_Change"][k])))
#             g.add((cnv_id_uri, rl.URIRef(impo + "has_tumor_sequence_allele_1"), rl.Literal(gene_mutations["Tumor_Seq_Allele1"][k])))
#             g.add((cnv_id_uri, rl.URIRef(impo + "has_tumor_sequence_allele_2"), rl.Literal(gene_mutations["Tumor_Seq_Allele2"][k])))
#             g.add((cnv_id_uri, rl.URIRef(impo + "has_start_position"), rl.Literal(gene_mutations["Start_position"][k])))
#             g.add((cnv_id_uri, rl.URIRef(impo + "has_end_position"), rl.Literal(gene_mutations["End_position"][k])))
#             g.add((cnv_id_uri, rl.URIRef(impo + "has_tumor_reference_count"), rl.Literal(gene_mutations["Tumor_ref_count"][k])))

# g.serialize(destination="transcriptomics_genomic_region_2.owl", format='xml') # generates owl file
# # del g
# # g = rl.Graph()
# # g.parse("transcriptomics_impo_7.owl")
# print("  ------  over 7")


##########################################################
#
#
#
#
#
################################################## BIOLOGICAL PATHWAY INSTANTIAION AND SAMPLE - BIOLOGICAL PATHWAY ASSIGNMENT

# for k in range(len(bio_pathway["Sample_ID"])):
#     sample_uri = rl.URIRef(placeholder +  bio_pathway["Sample_ID"][k].replace(" ","_"))
#     pathway_uri = rl.URIRef(placeholder + bio_pathway["biological_pathways"][k])
#     g.add((pathway_uri, RDF.type, rl.URIRef(impo + "biological_pathway"))) # adds as individual of biological pathway
#     sample_pathway_assignment = rl.URIRef(placeholder + "sample_biopath_id_" + bio_pathway["Sample_ID"][k] + "_" + bio_pathway["biological_pathways"][k])
#     g.add((sample_pathway_assignment, RDF.type, rl.URIRef(impo + "sample-biological_pathway_assignment")))
#     g.add((sample_pathway_assignment, rl.URIRef(impo + "has_source"), sample_uri))
#     g.add((sample_pathway_assignment, rl.URIRef(impo + "has_target"), pathway_uri))
#     g.add((sample_pathway_assignment, rl.URIRef(impo + "has_score"), rl.Literal(bio_pathway["biological_pathways_score"][k])))

# g.serialize(destination="transcriptomics_biopathway.owl", format='xml') # generates owl file
# # del g
# # g = rl.Graph()
# # g.parse("transcriptomics_impo_8.owl")
# print("  ------  over 8")


##########################################################
#
#
#
#
#
################################################## THERAPY

# ## on SUBJID file to make therapy and drug
# drugs = []
# for k in range(len(subjid["Arm"])):
#     drugs.append(subjid["Arm"][k])
# drugs = list(set(drugs)) # gets different kinds of drugs
# for i in range(len(drugs)):
#     drug_uri = rl.URIRef(placeholder + str(drugs[i])) # make drug individual uri
#     g.add((drug_uri, RDF.type, rl.URIRef(impo + "drug"))) # adds as individual of class drug
# for k in range(len(subjid["SUBJID"])):
#     therapy_uri = rl.URIRef(placeholder + "therapy_" + subjid["SUBJID"][k]) # makes therapy individual uri based on id of patient that received it
#     g.add((therapy_uri, RDF.type, rl.URIRef(impo + "therapy"))) # adds as individual of class therapy
#     g.add((rl.URIRef(placeholder + subjid["SUBJID"][k]), rl.URIRef(impo + "has_therapy"), therapy_uri)) # adds corresponding therapy to each patient
#     drug_uri = rl.URIRef(placeholder + str(subjid["Arm"][k]))
#     g.add((drug_uri, RDF.type, rl.URIRef(impo + "drug"))) # adds as individual of class drug
#     g.add((therapy_uri, rl.URIRef(impo + "has_drug"), drug_uri)) # adds drug used in therapy
#     g.add((therapy_uri, rl.URIRef(impo + "has_days_from_tumorsample_collection_and_start_of_trial_therapy"), rl.Literal(subjid["Days_from_TumorSample_Collection_and_Start_of_Trial_Therapy"][k]))) # adds Days_from_TumorSample_Collection_and_Start_of_Trial_Therapy
#     g.add((therapy_uri, rl.URIRef(impo + "has_benefit"), rl.Literal(subjid["Benefit"][k]))) # adds Benefit
#     g.add((therapy_uri, rl.URIRef(impo + "has_ORR"), rl.Literal(subjid["ORR"][k]))) # adds ORR
#     g.add((therapy_uri, rl.URIRef(impo + "has_PFS"), rl.Literal(subjid["PFS"][k]))) # adds PFS
#     g.add((therapy_uri, rl.URIRef(impo + "has_PFS_CNSR"), rl.Literal(subjid["PFS_CNSR"][k]))) # adds PFS_CNSR
#     g.add((therapy_uri, rl.URIRef(impo + "has_OS"), rl.Literal(subjid["OS"][k]))) # adds OS
#     g.add((therapy_uri, rl.URIRef(impo + "has_OS_CNSR"), rl.Literal(subjid["OS_CNSR"][k]))) # adds OS_CNSR

# g.serialize(destination="transcriptomics_therapy.owl", format='xml') # generates owl file
# # del g
# # g = rl.Graph()
# # g.parse("transcriptomics_impo_9.owl")
# print("  ------  over 9")



# ### INSTANTIATE CELL DECONVOLUTION AND CELL QUANTIFICATION CLASSES (CELL_DECONVOLUTION FILE)
# for k in range(len(cell_deconv["Sample_ID"])): # link sample to cell deconvolution assay id
#     for j in cell_deconv.columns:
#         cell_deconv_uri = rl.URIRef(placeholder + "assay_" + cell_deconv["Sample_ID"][k].replace(" ","_") + "_cell_deconv_" + j)
#         g.add((cell_deconv_uri, RDF.type, rl.URIRef(impo + "cell_deconvolution"))) # adds as individual of class cell deconv
#         cell_quantif_uri = rl.URIRef(placeholder + "cell_quantif_" + cell_deconv["Sample_ID"][k].replace(" ","_") + "_" + j)
#         g.add((cell_quantif_uri, RDF.type, rl.URIRef(impo + "cell_quantification"))) # adds as individual of class cell quantification
#         g.add((cell_quantif_uri, rl.URIRef(impo + "has_id"), rl.Literal(cell_deconv["Sample_ID"][k].replace(" ","_") + "_" + j)))
#         sample_uri = rl.URIRef(placeholder +  j.replace(" ","_"))
#         g.add((cell_quantif_uri, rl.URIRef(impo + "has_type"), sample_uri))
#         g.add((cell_quantif_uri, rl.URIRef(impo + "has_quantity"), rl.Literal(cell_deconv[j][k])))
#         g.add((cell_quantif_uri, rl.URIRef(impo + "identifies"), cell_deconv_uri)) # links cell quantification to cell deconvolution

# g.serialize(destination="transcriptomics_impo_cell_deconv.owl", format='xml') # generates owl file
# # del g
# # g = rl.Graph()
# # g.parse("transcriptomics_impo_9.owl")
# print("  ------  over 9")



subj_to_samp = {}
for k in range(len(meta["SUBJID"])):
    if subj_to_samp.get(meta["SUBJID"][k]) == None: subj_to_samp[meta["SUBJID"][k]] = []
    subj_to_samp[meta["SUBJID"][k]].append(meta["Sample_ID"][k])
    if subj_to_samp.get(gene_mutations["SUBJID"][k]) == None: subj_to_samp[gene_mutations["SUBJID"][k]] = []
    if gene_mutations["Sample_ID"][k] not in subj_to_samp[gene_mutations["SUBJID"][k]]:
        subj_to_samp[gene_mutations["SUBJID"][k]].append(gene_mutations["Sample_ID"][k])

# for k in range(len(tumor_feat["SUBJID"])):
#     if subj_to_samp.get(tumor_feat["SUBJID"][k]) != None:
#         for j in subj_to_samp[tumor_feat["SUBJID"][k]]:
#             sample_uri = rl.URIRef(placeholder + j.replace(" ","_"))
#             tissue_uri = rl.URIRef(placeholder + "tumor_" + j.replace(" ","_"))
#             g.add((tissue_uri, RDF.type, rl.URIRef(impo + "tumor_tissue"))) # adds as individual of class tumor tissue
#             g.add((sample_uri, rl.URIRef(impo + "has_source"), tissue_uri))
#             g.add((tissue_uri, rl.URIRef(impo + "has_ploidy"), rl.Literal(tumor_feat["Ploidy"][k])))
#             g.add((tissue_uri, rl.URIRef(impo + "has_purity"), rl.Literal(tumor_feat["Purity"][k])))
#             g.add((tissue_uri, rl.URIRef(impo + "has_immunophenotype"), rl.Literal(tumor_feat["ImmunoPhenotype"][k])))
#             g.add((tissue_uri, rl.URIRef(impo + "is_metastasis"), rl.Literal(tumor_feat["Tumor_Sample_Primary_or_Metastasis"][k])))
#             cancer_uri = rl.URIRef(placeholder + str(tumor_feat["Site_of_Metastasis"][k]).replace(" ","_"))
#             g.add((cancer_uri, RDF.type, rl.URIRef(impo + "cancer")))
#             g.add((tissue_uri, rl.URIRef(impo + "has_type"), rl.Literal(tumor_feat["Site_of_Metastasis"][k])))
#             g.add((tissue_uri, rl.URIRef(impo + "part_of"), cancer_uri))

# g.serialize(destination="transcriptomics_tumor_tissue.owl", format='xml') # generates owl file
# # del g
# # g = rl.Graph()
# # g.parse("transcriptomics_impo_10.owl")
# print("  ------  over 10")



for k in range(len(gene_mutations["SUBJID"])):  ### INSTANTIATE PROTEIN MUTATION AND ADD PROTEIN_CHANGE PROPERTY from GENE_MUTATIONS FILE
    try:
        sample_uri = rl.URIRef(placeholder + gene_mutations["Sample_ID"][k].replace(" ","_"))
        protein_mutation_uri = rl.URIRef(placeholder + gene_mutations["Sample_ID"][k].replace(" ","_") + "_" + str(gene_mutations["Protein_Change"][k]).replace(">","-"))
    except: # beacause gene_mutations file doesn't have all complete sample_id rows
        for j in subj_to_samp[gene_mutations["SUBJID"][k]]:
            sample_uri = rl.URIRef(placeholder + j.replace(" ","_"))
            protein_mutation_uri = rl.URIRef(placeholder + j.replace(" ","_") + "_" + str(gene_mutations["Protein_Change"][k]).replace(">","-"))
    g.add((protein_mutation_uri, RDF.type, rl.URIRef(impo + "protein_mutation"))) # adds as individual of class protein_mutation
    
    try:
        protein_uri = rl.URIRef(placeholder +  "gene_product_" + name_to_id[gene_mutations["gene_name"][k]][0])
        for m in subj_to_samp[gene_mutations["SUBJID"][k]]:
            genomic_mutation_uri = rl.URIRef(placeholder + "genomic_mutation_" + m)
            g.add((protein_mutation_uri, rl.URIRef(impo + "has_source"), genomic_mutation_uri))
            g.add((genomic_mutation_uri, rl.URIRef(impo + "has_source"), rl.URIRef(placeholder + name_to_id[gene_mutations["gene_name"][k]][0]))) ## ADD LINK BETWEEN GENE AND GENOMIC MUTATION??
        g.add((protein_mutation_uri, rl.URIRef(impo + "occurs_in"), protein_uri))
        g.add((protein_mutation_uri, rl.URIRef(impo + "has_protein_change"), rl.Literal(str(gene_mutations["Protein_Change"][k]))))
    except: pass

g.serialize(destination="transcriptomics_protein_mutation.owl", format='xml') # generates owl file

print("  ------  over -- THE END -- ")

