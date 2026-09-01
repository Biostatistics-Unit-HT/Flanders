# This script create an anndata starting from study and credisble_set parquet files from Open Target. 
# This anndata can later be used by the Flanders pipeline for colocalization 

import polars as pl
import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse import csr_matrix
import pyarrow as pa
import math
import anndata as ad
import argparse

#Although I give all the option here, some project have only certain type of data. For example FINNGEN_R12 has only gwas and SuSie
parser = argparse.ArgumentParser(description="Program to generate credible sets from Open Target in Anndata format for Flanders colocalization")
parser.add_argument("--project", default=None, type=str, help="Name of the project to get data from in OTG (example: FINNGEN_R12, GTEx, GCST)")
parser.add_argument("--type_sumstat", default='gwas', type=str, help="Either a gwas, eqtl or pqtl")
parser.add_argument("--type_finemap", default='SuSie', type=str, help="Either SuSie or SuSiE-inf")
parser.add_argument("--study", default=None, type=str, help="Parquet file of the study from OT")
parser.add_argument("--credset", default=None, type=str, help="Parquet file of the credible sets from OT")
parser.add_argument("--ref_panel", default=None, type=str, help="A csv panel to use to impute missing SNPs from the credible set. Mandatory columns are (chr,pos,a0,a1) wuth a0 and a1 in any particula order")
parser.add_argument("--out", default=None, type=str, help="Name of the output anndata file")
args = parser.parse_args()


class OT_create_anndata:
    def __init__(
        self, study:str, 
        project:str, 
        type_study:str, 
        credset:str, 
        ref_panel:str, 
        type_finemap:str):

        self.type_study = type_study
        self.study = study
        self.credset = credset
        self.project = project
        self.imputation_panel = ref_panel
        self.type_finemap = type_finemap
        self.log10 = math.log(10)
        self.log_prior_odds = math.log(0.01 / 0.99)

    def _read_study(self):
        self.study_pl = pl.read_parquet(self.study)
        self.study_pl = self.study_pl.filter((pl.col("hasSumstats") == True) & (pl.col("projectId") == self.project) & (pl.col('studyType') == self.type_study))
        self.study_pl = self.study_pl.with_columns(
            pl.col("ldPopulationStructure")
            .map_elements(
            lambda lst: ",".join(item["ldPopulation"] for item in lst) if lst is not None else None,
            return_dtype=pl.String  # avoid the warning
            )
            .alias("ldPopulationStructure")
            )
        if self.type_study != "eqtl" or self.type_study != "pqtl":
            self.study_pl = self.study_pl.with_columns(
                pl.col("traitFromSource").alias("biosampleId")
                )
        if self.project=='GTEx':
            self.study_pl = self.study_pl.filter(
                pl.col("studyId").str.contains("^gtex_ge_", literal=False)
            )
        if self.project=='UKB_PPP_EUR':
            self.study_pl = self.study_pl.filter(pl.col("ldPopulationStructure")=="nfe")
        return self.study_pl

    def _read_credset(self):
            print(self.type_study)
            credset_pl = pl.read_parquet(self.credset)
            if self.type_study == "eqtl" or self.type_study == "pqtl":
                self.credset_query_pl = credset_pl.filter((pl.col("isTransQtl") == False) & (pl.col('finemappingMethod') == self.type_finemap) & (pl.col("studyId").is_in(self.study_pl["studyId"].to_list())))
            else:
                self.credset_query_pl = credset_pl.filter( (pl.col('finemappingMethod') == self.type_finemap) & (pl.col("studyId").is_in(self.study_pl["studyId"].to_list())))
                print('ok')
            self.credset_query_pl = self.credset_query_pl.filter(pl.col("locus").list.eval(pl.element().struct.field("is99CredibleSet")).list.any())
            if len(self.credset_query_pl)==0:
                print('Empty credible sets, please check if the parameter are right (type_sumstat, project, type_finemap)')
                exit()
            

    def _reshape_credset(self):
        study_credset = self.credset_query_pl.join(self.study_pl, on = "studyId", how = "inner")
        if "biosampleId" not in study_credset.columns:
            study_credset = study_credset.with_columns(pl.lit(None).alias("biosampleId"))
        var = study_credset.explode("locus").with_columns(pl.col("locus").struct.unnest())
        var = (var
                .group_by("studyLocusId")
                .agg([
                    pl.col("logBF").sum().alias("sumlogABF")
                ])
        )
        study_credset = var.join(study_credset, on = "studyLocusId")
        study_credset = study_credset.with_columns(
                pl.when((pl.lit(self.type_study != "eqtl")) | (pl.lit(self.type_study != "pqtl"))).then(
                pl.col("studyId").alias("phenotype_id")
                ).otherwise(
                    pl.concat_str(
                    [
                pl.col("biosampleId"),
                pl.lit(":"),
                pl.col("geneId")
                    ]
                ).alias("phenotype_id")
                ))
        study_credset = study_credset.with_columns(
                    pl.concat_str([
                        pl.col("projectId"),
                        pl.lit(":"),
                        pl.col("studyType")
                    ]).alias("projectId")  
                )
        study_credset = study_credset.with_columns([pl.col('variantId').alias('snp'), pl.col('chromosome').alias('chr')])
        study_credset = study_credset.with_columns(
            pl.col("snp").str.split_exact("_", 4)
            .struct.rename_fields(["chr_snp", "pos_snp", "a0", "a1"])
            .alias("fields")
        ).unnest("fields").drop("snp").with_columns([
            pl.min_horizontal("a0", "a1").alias("a_allele1"),
            pl.max_horizontal("a0", "a1").alias("a_allele2")
        ]).with_columns(
            pl.concat_str([
                pl.lit('chr'),
                pl.col("chr_snp").cast(pl.Utf8).str.replace("X", "23"),
                pl.lit(":"),
                pl.col("pos_snp"),
                pl.lit(":"),
                pl.col("a_allele1"),
                pl.lit(":"),
                pl.col("a_allele2"),
            ]).alias("snp")
        ).drop(["chr_snp","pos_snp","a_allele2","a_allele1","a0","a1"])
        study_credset = study_credset.with_columns(pl.col("chr").cast(str).replace("X", "23").alias("chr"))
        self.study_credset = study_credset.with_columns(
            pl.concat_str([
                pl.lit('chr'),
                pl.col("chr"),
                pl.lit("::"),
                pl.col("projectId"),
                pl.lit("::"),
                pl.col("phenotype_id"),
                pl.lit("::"),
                pl.col("snp")
            ]).alias("csname"),
            pl.when((pl.lit(self.type_study == "eqtl")) | (pl.lit(self.type_study == "pqtl"))).then(
            (pl.col("region").str.extract(r":(-?\d+)", 1).cast(pl.Int64) - 500000).clip(lower_bound = 0)
            ).otherwise(
            pl.col("region").str.extract(r":(-?\d+)", 1).cast(pl.Int64)
            ).alias("start"),
            pl.when(pl.lit(self.type_study == "eqtl") | (pl.lit(self.type_study == "pqtl"))).then(
            pl.col("region").str.extract(r"-(-?\d+)$", 1).cast(pl.Int64) + 500000
            ).otherwise(
            pl.col("region").str.extract(r"-(-?\d+)$", 1).cast(pl.Int64)
            ).alias("end"),
            (pl.col("sumlogABF") * (self.log10 + self.log_prior_odds - (pl.col("locus").list.len().cast(pl.Float64).log()))).alias("min_res_labf"),
            pl.col("sumlogABF").alias("cs_log10bf"),
            pl.col('snp').alias("topsnp")
            ).drop("chr")

    def _explode_locus_snps(self):
        if "traitFromSource" not in self.study_credset.columns:
            self.study_credset = self.study_credset.with_columns(pl.lit(None).alias("traitFromSource"))
        var = self.study_credset.explode("locus").with_columns(pl.col("locus").struct.unnest())
        var= var.with_columns(
                pl.col("variantId").str.split_exact("_", 4)
                .struct.rename_fields(["chr", "pos", "a0", "a1"])
                .alias("fields")
            ).unnest("fields").drop("variantId").with_columns([
            pl.min_horizontal("a0", "a1").alias("a_allele1"),
            pl.max_horizontal("a0", "a1").alias("a_allele2")
            ]).with_columns(
                pl.concat_str([
                    pl.lit("chr"),
                    pl.col("chr").cast(pl.Utf8).str.replace("X", "23"),  # Note: you have "chr" twice - is this intentional?
                    pl.lit(":"),
                    pl.col("pos"),
                pl.lit(":"),
                pl.col("a_allele1"),
                pl.lit(":"),
                pl.col("a_allele2")
            ]).alias("snp"),
            pl.concat_str([
                pl.lit("chr"),
                pl.col("chr")
            ]).alias("chrname"),
            pl.col('standardError').alias('se')
        ).drop("chr","standardError")
        var = var.with_columns(pl.col("chrname").cast(str).replace("chrX", "chr23").alias("chr"))
        var = var.with_columns(pl.col("pos").cast(pl.Int64)).drop("chrname")
        self.var = var.with_columns([
            pl.when(
            (pl.col("a0") == pl.col("a_allele2")) & (pl.col("a1") == pl.col("a_allele1"))
            ).then(
                -pl.col("beta")
            ).otherwise(
                pl.col("beta")
            ).alias("beta")]).select(
                "csname",
                'projectId',
                'start',
                'end',
                'phenotype_id',
                'topsnp',
                'min_res_labf',
                'cs_log10bf',
                'logBF',
                'beta',
                'se',
                'snp',
                'chr',
                'pos',
                'a_allele1',
                'a_allele2',
                'nSamples'
            )
    
    def _merge_with_imputation_refined(self):
        imputation_panel = pl.read_csv(self.imputation_panel)
        ref = imputation_panel.with_columns([
            pl.min_horizontal("a0", "a1").alias("a_allele1"),
            pl.max_horizontal("a0", "a1").alias("a_allele2")
        ]).with_columns(
            pl.col("chr").cast(pl.Utf8).str.replace("chrX", "chr23"),
            pl.concat_str([
            pl.col("chr").cast(pl.Utf8).str.replace("chrX", "chr23"),
            pl.lit(":"),
            pl.col("pos"),
            pl.lit(":"),
            pl.col("a_allele1"),
            pl.lit(":"),
            pl.col("a_allele2")
        ]).alias("snp")
        )
    
    # Perform full outer join
        ot_merge = ref.join(
        self.var,
        on=['chr', 'pos', 'a_allele1', 'a_allele2'],
        how='full'
        )
    
    # Coalesce the columns
        ot_merge = ot_merge.with_columns([
        pl.coalesce(["chr", "chr_right"]).alias("chr"),
        pl.coalesce(["pos", "pos_right"]).alias("pos"),
        pl.coalesce(["a_allele1", "a_allele1_right"]).alias("a_allele1"),
        pl.coalesce(["a_allele2", "a_allele2_right"]).alias("a_allele2"),
        pl.coalesce(["snp", "snp_right"]).alias("snp")
        ]).drop(["chr_right", "pos_right", "a_allele1_right", "a_allele2_right", "snp_right"])
    
        # Sort by chromosome and position
        sorted_df = ot_merge.sort(["chr", "pos"])
    
        # Create a helper column to identify credible set boundaries
        # We'll mark the start of each credible set region
        result_with_filled_nulls = sorted_df.with_columns([
        # Create boundary markers for credible sets
        pl.when(
            pl.col("csname").is_not_null() & 
            pl.col("start").is_not_null() & 
            pl.col("end").is_not_null()
        ).then(1).otherwise(0).alias("is_cs_boundary"),
        
        # Forward fill the credible set information within chromosomal regions
        pl.col("csname").forward_fill().over("chr"),
        pl.col("start").forward_fill().over("chr"),
        pl.col("end").forward_fill().over("chr"),
        pl.col("projectId").forward_fill().over("chr"),
        pl.col("phenotype_id").forward_fill().over("chr"),
        pl.col("topsnp").forward_fill().over("chr"),
        pl.col("min_res_labf").forward_fill().over("chr"),
        pl.col("cs_log10bf").forward_fill().over("chr"),
        ]).with_columns([
        # Only keep the filled values if the SNP falls within the credible set region
        pl.when(
            (pl.col("pos") >= pl.col("start")) & 
            (pl.col("pos") <= pl.col("end"))
        ).then(pl.col("csname")).otherwise(None).alias("csname"),
        
        pl.when(
            (pl.col("pos") >= pl.col("start")) & 
            (pl.col("pos") <= pl.col("end"))
        ).then(pl.col("projectId")).otherwise(None).alias("projectId"),
        
        pl.when(
            (pl.col("pos") >= pl.col("start")) & 
            (pl.col("pos") <= pl.col("end"))
        ).then(pl.col("phenotype_id")).otherwise(None).alias("phenotype_id"),
        
        pl.when(
            (pl.col("pos") >= pl.col("start")) & 
            (pl.col("pos") <= pl.col("end"))
        ).then(pl.col("topsnp")).otherwise(None).alias("topsnp"),
        
        pl.when(
            (pl.col("pos") >= pl.col("start")) & 
            (pl.col("pos") <= pl.col("end"))
        ).then(pl.col("min_res_labf")).otherwise(None).alias("min_res_labf"),
        
        pl.when(
            (pl.col("pos") >= pl.col("start")) & 
            (pl.col("pos") <= pl.col("end"))
        ).then(pl.col("cs_log10bf")).otherwise(None).alias("cs_log10bf"),
        
        pl.when(
            (pl.col("pos") >= pl.col("start")) & 
            (pl.col("pos") <= pl.col("end"))
        ).then(pl.col("start")).otherwise(None).alias("start"),
        
        pl.when(
            (pl.col("pos") >= pl.col("start")) & 
            (pl.col("pos") <= pl.col("end"))
        ).then(pl.col("end")).otherwise(None).alias("end"),
        ]).with_columns([
        # Fill summary statistics with 0
        pl.col("logBF").fill_null(0),
        pl.col("beta").fill_null(0),
        pl.col("se").fill_null(0),
        
        # Ensure SNP format is consistent
        pl.col("snp").fill_null(
            pl.concat_str([
                pl.col("chr"),
                pl.lit(":"),
                pl.col("pos"),
                pl.lit(":"),
                pl.col("a_allele1"),
                pl.lit(":"),
                pl.col("a_allele2"),
            ])
        )
        ]).select([
        "chr", 
        pl.col("pos").cast(pl.Utf8), 
        "csname",
        "snp",    
        'start',
        'end',
        'projectId',
        'phenotype_id',
        'topsnp',
        'min_res_labf',
        'cs_log10bf',
        'logBF',
        'beta',
        'se'
        ])
    
        self.result_with_filled_nulls = result_with_filled_nulls
    
    
    def _create_anndata_refined(self):
        susie_credset_pd = self.result_with_filled_nulls.to_pandas()
    
        # Only use rows with valid credible sets
        valid_data = susie_credset_pd[susie_credset_pd["csname"].notna()].copy()
        valid_data["csname"] = valid_data["csname"].astype(str)
    
        # Factorize using only valid data to ensure consistent lengths
        row_idx, row_labels = valid_data["csname"].factorize()
        col_idx, col_labels = valid_data["snp"].factorize()
    
        # Get values from the same valid_data to ensure consistent lengths
        logbf = valid_data["logBF"].fillna(0).to_numpy()
        beta = valid_data["beta"].fillna(0).to_numpy()
        se = valid_data["se"].fillna(0).to_numpy()
    
        # Create sparse matrices
        sparse_matrix_logbf = coo_matrix((logbf, (row_idx, col_idx)))
        sparse_matrix_beta = coo_matrix((beta, (row_idx, col_idx)))
        sparse_matrix_se = coo_matrix((se, (row_idx, col_idx)))
    
        adata = ad.AnnData(sparse_matrix_logbf)
        adata.obs["cs_name"] = row_labels
        adata.var["snp"] = col_labels
        adata.X = csr_matrix(adata.X)
        adata.layers["beta"] = csr_matrix(sparse_matrix_beta)
        adata.layers["se"] = csr_matrix(sparse_matrix_se)
    
        # Create mappings for obs metadata
        start_mapping = dict(zip(valid_data['csname'], valid_data['start']))
        end_mapping = dict(zip(valid_data['csname'], valid_data['end']))
        topsnp_mapping = dict(zip(valid_data['csname'], valid_data['topsnp']))
        credset_abf = dict(zip(valid_data['csname'], valid_data['cs_log10bf']))
        project_mapping = dict(zip(valid_data['csname'], valid_data['projectId']))
        min_credset_abf = dict(zip(valid_data['csname'], valid_data['min_res_labf']))
        study_mapping = dict(zip(valid_data['csname'], valid_data['phenotype_id']))
    
        # Assign obs metadata
        start_labels = [start_mapping.get(label) for label in row_labels]
        end_labels = [end_mapping.get(label) for label in row_labels]
        topsnp_labels = [topsnp_mapping.get(label) for label in row_labels]
        credset_abf_labels = [credset_abf.get(label) for label in row_labels]
        phenotype_labels = [study_mapping.get(label) for label in row_labels]
        min_credset_abf_labels = [min_credset_abf.get(label) for label in row_labels]
        project_labels = [project_mapping.get(label) for label in row_labels]
    
        adata.obs["start"] = start_labels
        adata.obs["end"] = end_labels
        adata.obs["topsnp"] = topsnp_labels
        adata.obs["study_id"] = project_labels
        adata.obs["min_res_labf"] = min_credset_abf_labels
        adata.obs["logsum(logABF)"] = credset_abf_labels
        adata.obs['start'] = adata.obs['start'].apply(lambda x: 0 if x < 0 else x)
        adata.obs["phenotype_id"] = phenotype_labels
        adata.obs[['chr']] = adata.obs['topsnp'].str.split(':', expand=True)[[0]]
        adata.obs.loc[:, "coverage"] = 0.99
    
        # Assign var metadata
        adata.var[['chr', 'pos']] = adata.var['snp'].str.split(':', expand=True)[[0, 1]]
    
        # Set indices
        adata.var.index = adata.var['snp']
        adata.obs.index = adata.obs['cs_name']    
        return adata
    
if __name__=='__main__':
    ot_obj = OT_create_anndata(study = args.study, credset = args.credset, ref_panel = args.ref_panel, type_finemap=args.type_finemap, project=args.project, type_study = args.type_sumstat)
    ot_obj._read_study()
    ot_obj._read_credset()
    ot_obj._reshape_credset()
    ot_obj._explode_locus_snps()
    ot_obj._merge_with_imputation_refined()
    adata = ot_obj._create_anndata()
    adata.write(f'{args.out}.h5ad', compression="gzip")

