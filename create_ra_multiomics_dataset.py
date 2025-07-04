#!/usr/bin/env python3
"""
Generate synthetic multiomics dataset for rheumatoid arthritis analysis
Based on the methodology from Alsaedi et al. (2024)
"""

import pandas as pd
import numpy as np
import random
from datetime import datetime, timedelta

# Set random seed for reproducibility
np.random.seed(42)
random.seed(42)

print("Generating Rheumatoid Arthritis Multiomics Dataset...")
print("=" * 60)

# Define RA-associated genes from literature
ra_genes = [
    # HLA genes (major histocompatibility complex)
    "HLA-DRB1", "HLA-DQA1", "HLA-DQB1", "HLA-DPA1", "HLA-DPB1",
    
    # Cytokine and immune response genes
    "TNF", "IL1B", "IL6", "IL10", "IL17A", "IL17F", "IL23R", "IL2RA",
    "IFNG", "IL4", "IL13", "IL12B", "IL18", "IL1RN", "TNFAIP3",
    
    # T cell and B cell related genes
    "CD40", "CD80", "CD86", "CTLA4", "PDCD1", "CD28", "ICOS", "CD19",
    "CD20", "CD22", "BAFF", "APRIL", "TACI", "BCMA",
    
    # Transcription factors and signaling
    "STAT1", "STAT3", "STAT4", "STAT6", "IRF4", "IRF5", "IRF7", "IRF8",
    "NFKB1", "NFKB2", "REL", "RELA", "RELB", "FOXP3", "TBX21", "GATA3",
    
    # Complement and coagulation
    "C1QA", "C1QB", "C1QC", "C3", "C4A", "C4B", "CFB", "CFD", "CFH",
    "F2", "F3", "F5", "F7", "F8", "F10", "SERPINE1", "PLAT", "PLG",
    
    # Matrix metalloproteinases and tissue remodeling
    "MMP1", "MMP3", "MMP9", "MMP13", "TIMP1", "TIMP2", "COL1A1", "COL2A1",
    
    # Chemokines and receptors
    "CCL2", "CCL3", "CCL4", "CCL5", "CCL20", "CXCL8", "CXCL10", "CXCL12",
    "CCR1", "CCR2", "CCR5", "CXCR3", "CXCR4", "CXCR5",
    
    # Apoptosis and cell cycle
    "BCL2", "BAX", "FAS", "FASLG", "CASP3", "CASP8", "TP53", "MDM2",
    
    # Bone and cartilage metabolism
    "RANKL", "OPG", "RANK", "RUNX2", "SOX9", "ACAN", "COL2A1", "COMP",
    
    # Additional RA-associated genes
    "PTPN22", "PADI4", "FCGR3A", "TRAF1", "C5", "BLK", "BANK1", "TNIP1",
    "PRKCQ", "KIF5A", "TAGAP", "RASGRP1", "AFF3", "CD247", "SPRED2"
]

# Pathway categories for RA
pathway_categories = {
    'Immune_Response': [
        'HLA-DRB1', 'HLA-DQA1', 'STAT1', 'STAT3', 'STAT4', 'IRF4', 'IRF5',
        'NFKB1', 'RELA', 'CD40', 'CD80', 'CD86', 'CTLA4', 'IL2RA'
    ],
    'Inflammatory_Response': [
        'TNF', 'IL1B', 'IL6', 'IL17A', 'IL23R', 'IFNG', 'TNFAIP3',
        'CCL2', 'CCL3', 'CXCL8', 'CXCL10', 'IL1RN'
    ],
    'T_Cell_Activation': [
        'CD28', 'ICOS', 'FOXP3', 'TBX21', 'GATA3', 'IL4', 'IL13',
        'PDCD1', 'CD247', 'PRKCQ', 'RASGRP1'
    ],
    'B_Cell_Function': [
        'CD19', 'CD20', 'CD22', 'BAFF', 'APRIL', 'TACI', 'BCMA',
        'BLK', 'BANK1', 'AFF3'
    ],
    'Bone_Cartilage_Metabolism': [
        'RANKL', 'OPG', 'RANK', 'RUNX2', 'SOX9', 'ACAN', 'COL1A1',
        'COL2A1', 'COMP', 'MMP1', 'MMP3', 'MMP9', 'MMP13', 'TIMP1'
    ],
    'Complement_Coagulation': [
        'C1QA', 'C3', 'C4A', 'C5', 'CFB', 'CFH', 'F2', 'F3', 'F8',
        'SERPINE1', 'PLAT', 'PLG'
    ],
    'Apoptosis_Cell_Cycle': [
        'BCL2', 'BAX', 'FAS', 'FASLG', 'CASP3', 'CASP8', 'TP53', 'MDM2'
    ]
}

# Generate genetic variants dataset
print("1. Generating genetic variants dataset...")
n_variants = 150  # More variants than the COVID study
genetic_variants = []

for i, gene in enumerate(ra_genes):
    # Generate 1-3 variants per gene
    n_variants_per_gene = np.random.choice([1, 2, 3], p=[0.5, 0.3, 0.2])
    
    for j in range(n_variants_per_gene):
        variant = {
            'variant_id': f'rs{np.random.randint(1000000, 99999999)}',
            'chromosome': f'chr{np.random.randint(1, 23)}',
            'position': np.random.randint(1000000, 200000000),
            'gene_symbol': gene,
            'ref_allele': np.random.choice(['A', 'T', 'G', 'C']),
            'alt_allele': np.random.choice(['A', 'T', 'G', 'C']),
            'p_value': np.random.uniform(1e-12, 5e-6),  # Stronger associations
            'odds_ratio': np.random.uniform(1.2, 4.5),
            'allele_frequency': np.random.uniform(0.05, 0.45),
            'functional_consequence': np.random.choice([
                'missense_variant', 'synonymous_variant', 'intron_variant',
                'upstream_gene_variant', 'downstream_gene_variant', 
                '3_prime_UTR_variant', '5_prime_UTR_variant', 'splice_site_variant'
            ]),
            'population': np.random.choice(['European', 'East_Asian', 'African', 'Mixed'], 
                                        p=[0.4, 0.3, 0.2, 0.1]),
            'study_type': np.random.choice(['GWAS', 'Candidate_gene', 'Exome_seq', 'WGS']),
            'sample_size': np.random.randint(1000, 50000),
            'effect_direction': np.random.choice(['risk', 'protective'], p=[0.8, 0.2])
        }
        genetic_variants.append(variant)

genetic_df = pd.DataFrame(genetic_variants)
print(f"   Generated {len(genetic_df)} genetic variants for {len(ra_genes)} genes")

# Generate gene expression dataset
print("2. Generating gene expression dataset...")
n_samples = 200  # 100 RA patients, 100 controls
conditions = ['RA_patient'] * 100 + ['healthy_control'] * 100
tissues = np.random.choice(['synovial_tissue', 'peripheral_blood', 'synovial_fluid'], 
                          n_samples, p=[0.4, 0.4, 0.2])

expression_data = []
for i in range(n_samples):
    sample_data = {
        'sample_id': f'SAMPLE_{i+1:03d}',
        'condition': conditions[i],
        'tissue_type': tissues[i],
        'age': np.random.randint(25, 80),
        'sex': np.random.choice(['Male', 'Female']),
        'disease_duration': np.random.randint(0, 20) if conditions[i] == 'RA_patient' else 0,
        'treatment_status': np.random.choice(['naive', 'DMARD', 'biologic', 'combination']) 
                           if conditions[i] == 'RA_patient' else 'none'
    }
    
    # Generate expression values for each gene
    for gene in ra_genes:
        # RA patients generally have altered expression
        if conditions[i] == 'RA_patient':
            # Some genes are upregulated, others downregulated
            if gene in ['TNF', 'IL1B', 'IL6', 'IL17A', 'MMP1', 'MMP3', 'RANKL']:
                # Upregulated in RA
                expression = np.random.normal(8.5, 1.2)  # Higher expression
            elif gene in ['IL10', 'FOXP3', 'OPG', 'TIMP1', 'IL1RN']:
                # Downregulated in RA
                expression = np.random.normal(5.5, 1.0)  # Lower expression
            else:
                # Moderately altered
                expression = np.random.normal(7.0, 1.5)
        else:
            # Healthy controls have more normal expression
            expression = np.random.normal(6.8, 1.0)
        
        sample_data[f'{gene}_expression'] = max(0, expression)  # Ensure non-negative
    
    expression_data.append(sample_data)

expression_df = pd.DataFrame(expression_data)
print(f"   Generated expression data for {n_samples} samples and {len(ra_genes)} genes")

# Generate protein abundance dataset
print("3. Generating protein abundance dataset...")
protein_data = []
for i in range(n_samples):
    sample_data = {
        'sample_id': f'SAMPLE_{i+1:03d}',
        'condition': conditions[i],
        'tissue_type': tissues[i]
    }
    
    # Generate protein abundance values
    for gene in ra_genes:
        # Protein levels generally correlate with mRNA but with some noise
        base_expression = expression_df.loc[i, f'{gene}_expression']
        
        # Add biological noise and post-translational effects
        protein_abundance = base_expression + np.random.normal(0, 0.5)
        
        # Some proteins have different regulation than mRNA
        if gene in ['TNF', 'IL1B', 'IL6'] and conditions[i] == 'RA_patient':
            protein_abundance += np.random.normal(1.0, 0.3)  # Enhanced protein production
        
        sample_data[f'{gene}_protein'] = max(0, protein_abundance)
    
    protein_data.append(sample_data)

protein_df = pd.DataFrame(protein_data)
print(f"   Generated protein data for {n_samples} samples and {len(ra_genes)} proteins")

# Generate metabolomics dataset
print("4. Generating metabolomics dataset...")
# Define metabolites associated with RA
ra_metabolites = [
    'Glucose', 'Lactate', 'Pyruvate', 'Citrate', 'Succinate', 'Fumarate',
    'Alanine', 'Glycine', 'Serine', 'Threonine', 'Valine', 'Leucine',
    'Isoleucine', 'Proline', 'Phenylalanine', 'Tyrosine', 'Tryptophan',
    'Histidine', 'Arginine', 'Lysine', 'Glutamate', 'Glutamine', 'Aspartate',
    'Asparagine', 'Cysteine', 'Methionine', 'Creatine', 'Creatinine',
    'Urea', 'Uric_acid', 'Cholesterol', 'Triglycerides', 'Palmitic_acid',
    'Stearic_acid', 'Oleic_acid', 'Linoleic_acid', 'Arachidonic_acid',
    'Prostaglandin_E2', 'Leukotriene_B4', 'Thromboxane_B2'
]

metabolomics_data = []
for i in range(n_samples):
    sample_data = {
        'sample_id': f'SAMPLE_{i+1:03d}',
        'condition': conditions[i],
        'tissue_type': tissues[i]
    }
    
    for metabolite in ra_metabolites:
        # RA patients have altered metabolite profiles
        if conditions[i] == 'RA_patient':
            if metabolite in ['Lactate', 'Glucose', 'Arachidonic_acid', 'Prostaglandin_E2']:
                # Elevated in RA
                concentration = np.random.lognormal(2.5, 0.8)
            elif metabolite in ['Glutamine', 'Arginine', 'Cholesterol']:
                # Decreased in RA
                concentration = np.random.lognormal(1.8, 0.6)
            else:
                concentration = np.random.lognormal(2.2, 0.7)
        else:
            concentration = np.random.lognormal(2.0, 0.5)
        
        sample_data[f'{metabolite}_concentration'] = concentration
    
    metabolomics_data.append(sample_data)

metabolomics_df = pd.DataFrame(metabolomics_data)
print(f"   Generated metabolomics data for {n_samples} samples and {len(ra_metabolites)} metabolites")

# Generate clinical data
print("5. Generating clinical data...")
clinical_data = []
for i in range(n_samples):
    if conditions[i] == 'RA_patient':
        clinical_record = {
            'sample_id': f'SAMPLE_{i+1:03d}',
            'condition': conditions[i],
            'age': expression_df.loc[i, 'age'],
            'sex': expression_df.loc[i, 'sex'],
            'disease_duration_years': expression_df.loc[i, 'disease_duration'],
            'RF_positive': np.random.choice([True, False], p=[0.7, 0.3]),
            'ACPA_positive': np.random.choice([True, False], p=[0.65, 0.35]),
            'DAS28_score': np.random.uniform(2.6, 7.8),  # Disease Activity Score
            'HAQ_score': np.random.uniform(0, 3.0),      # Health Assessment Questionnaire
            'ESR_mm_hr': np.random.uniform(15, 100),     # Erythrocyte Sedimentation Rate
            'CRP_mg_L': np.random.uniform(3, 50),        # C-Reactive Protein
            'joint_count_swollen': np.random.randint(0, 28),
            'joint_count_tender': np.random.randint(0, 28),
            'morning_stiffness_min': np.random.randint(30, 240),
            'treatment_response': np.random.choice(['good', 'moderate', 'poor'], p=[0.4, 0.4, 0.2])
        }
    else:
        clinical_record = {
            'sample_id': f'SAMPLE_{i+1:03d}',
            'condition': conditions[i],
            'age': expression_df.loc[i, 'age'],
            'sex': expression_df.loc[i, 'sex'],
            'disease_duration_years': 0,
            'RF_positive': False,
            'ACPA_positive': False,
            'DAS28_score': np.random.uniform(0, 2.5),
            'HAQ_score': np.random.uniform(0, 0.5),
            'ESR_mm_hr': np.random.uniform(2, 20),
            'CRP_mg_L': np.random.uniform(0.1, 3.0),
            'joint_count_swollen': 0,
            'joint_count_tender': 0,
            'morning_stiffness_min': np.random.randint(0, 30),
            'treatment_response': 'not_applicable'
        }
    
    clinical_data.append(clinical_record)

clinical_df = pd.DataFrame(clinical_data)
print(f"   Generated clinical data for {n_samples} samples")

# Generate pathway annotation data
print("6. Generating pathway annotation data...")
pathway_data = []
for pathway, genes in pathway_categories.items():
    for gene in genes:
        pathway_record = {
            'gene_symbol': gene,
            'pathway_name': pathway,
            'pathway_category': pathway.replace('_', ' '),
            'gene_function': f'Function related to {pathway.replace("_", " ").lower()}',
            'pathway_size': len(genes),
            'pathway_significance': np.random.uniform(1e-8, 1e-3)
        }
        pathway_data.append(pathway_record)

pathway_df = pd.DataFrame(pathway_data)
print(f"   Generated pathway annotations for {len(pathway_df)} gene-pathway associations")

# Save all datasets
print("\n7. Saving datasets...")
genetic_df.to_csv('/home/ubuntu/ra_genetic_variants.csv', index=False)
expression_df.to_csv('/home/ubuntu/ra_gene_expression.csv', index=False)
protein_df.to_csv('/home/ubuntu/ra_protein_abundance.csv', index=False)
metabolomics_df.to_csv('/home/ubuntu/ra_metabolomics.csv', index=False)
clinical_df.to_csv('/home/ubuntu/ra_clinical_data.csv', index=False)
pathway_df.to_csv('/home/ubuntu/ra_pathway_annotations.csv', index=False)

# Generate summary statistics
print("\n" + "=" * 60)
print("DATASET SUMMARY")
print("=" * 60)
print(f"Genetic Variants: {len(genetic_df)} variants across {genetic_df['gene_symbol'].nunique()} genes")
print(f"Gene Expression: {len(expression_df)} samples, {len(ra_genes)} genes")
print(f"Protein Abundance: {len(protein_df)} samples, {len(ra_genes)} proteins")
print(f"Metabolomics: {len(metabolomics_df)} samples, {len(ra_metabolites)} metabolites")
print(f"Clinical Data: {len(clinical_df)} samples")
print(f"Pathway Annotations: {len(pathway_df)} gene-pathway associations")
print(f"Pathways: {len(pathway_categories)} major biological pathways")

print(f"\nCondition Distribution:")
print(f"  RA Patients: {sum(clinical_df['condition'] == 'RA_patient')}")
print(f"  Healthy Controls: {sum(clinical_df['condition'] == 'healthy_control')}")

print(f"\nTissue Distribution:")
tissue_counts = clinical_df.groupby('condition')['sample_id'].count()
print(f"  Synovial Tissue: {sum(expression_df['tissue_type'] == 'synovial_tissue')}")
print(f"  Peripheral Blood: {sum(expression_df['tissue_type'] == 'peripheral_blood')}")
print(f"  Synovial Fluid: {sum(expression_df['tissue_type'] == 'synovial_fluid')}")

print("\nDatasets saved successfully!")
print("Ready for multiomics network analysis tutorial!")

