## list of species


See Methods from raw files download for the sources.

Species abreviation| Binomial name                    | Source
-------------------|----------------------------------|-----------------------
ACAPL              | Acanthaster planci               | NCBI GCF_001949145.1
ACRMI              | Acropora millepora               | NCBI GCF_004143615.1
ADIVA              | Adineta vaga                     | ENSEMBL Release-45
AMPQU              | Amphimedon queenslandica         | ENSEMBL Release-45
ANOGA              | Anopheles gambiae                | ENSEMBL Release-45
AURAU              | Aurelia aurita                   | David Gold Google Drive Aurelia.Genome_v1.2        
BRALA              | Branchiostoma lanceolatum        | ENSEMBL Release-45
CAEEL              | Caenorhabditis elegans           | ENSEMBL Release-45              
CALMI              | Callorhinchus milii              | NCBI GCF_000165045.1
CAPOW              | Capitella teleta                 | ENSEMBL Release-45
CAPTE              | Capsaspora owczarzaki            | ENSEMBL Release-45
CHEMY              | Chelonia mydas                   | NCBI GCF_000344595.1
CIOIN              | Ciona intestinalis               | NCBI GCF_000224145.3
CLYHE              | Clytia hemisphaerica             | MARIMBA
CRAGI              | Crassostrea gigas                | ENSEMBL Release-45 
DANRE              | Danio rerio                      | NCBI GCF_000002035.6
DAPPU              | Daphnia pulex                    | ENSEMBL Release-45
DROME              | Drosophila melanogaster          | ENSEMBL Release-45
EUPSC              | Euprymna scolopes                | Lachesis assembly (Schmidbaur et al. in prep.)
EXAPA              | Exaiptasia pallida               | NCBI GCF_001417965.1
GALGA              | Gallus gallus                    | NCBI GCF_000002315.6
HELRO              | Helobdella robusta               | ENSEMBL Release-45
HIPCO              | Hippocampus comes                | NCBI GCF_001891065.1
HOFMI              | Hofstenia miamia                 | Downloaded from http://srivastavalab.rc.fas.harvard.edu also available in ENSEMBL Release-45
HOIHO              | Hoilungia hongkongensis          | https://bitbucket.org/molpalmuc/hoilungia-genome/src/master/tracks/
HOMSA              | Homo sapiens                     | NCBI GCF_000001405.39         
HYDVU              | Hydra vulgaris                   | NHGRI hydra2.0
IXOSC              | Ixodes scapularis                | ENSEMBL Release-45
LATCH              | Latimeria chalumnae              | NCBI GCF_000225785.1
LEPOC              | Lepisosteus oculatus             | NCBI GCF_000242695.1      
LINAN              | Lingula anatina                  | ENSEMBL Release-45
LOTGI              | Lottia gigantea                  | ENSEMBL Release-45
MAYZE              | Maylandia zebra                  | NCBI GCF_000238955.4      
MIZYE              | Mizuhopecten yessoensis          | Wang et al. 2017
MNELE              | Mnemiopsis leidyi                | NHGRI ML2.2
MUSMU              | Mus musculus                     | NCBI GCF_000001635.26
NEMVE              | Nematostella vectensis           | ENSEMBL Release-45
PARTE              | Parasteatoda tepidariorum        | NCBI GCF_000365465.2 
PLEBA              | Pleurobrachia bachei             | NCBI GCA_000695325.1, Neurobase      
PTYFL              | Ptychodera flava                 | OIST pfl_public_ver1.0
SACKO              | Saccoglossus kowalevskii         | OIST Sackov3
SALRO              | Salpingoeca rosetta              | ENSEMBL Release-45
SCHME              | Schmidtea mediterranea           | Planmine
STRMA              | Strigamia maritima               | ENSEMBL Release-45
STRPU              | Strongylocentrotus purpuratus    | ENSEMBL Release-45
SYCCI              | Sycon ciliatum                   | DataDryad
TRIAD              | Trichoplax adhaerens             | ENSEMBL Release-45     
TRICA              | Tribolium castaneum              | ENSEMBL Release-45
XENTR              | Xenopus tropicalis               | NCBI GCF_000004195.3

## External scripts
orthoFinderToOrthogroup.pl: https://github.com/nijibabulu/metazoan_synteny/tree/master/scripts

## 1. Download of the raw data
### 1.1. Download peptide files
42/49 peptide databases directly downloaded using curl (or copied from proj, or cloned with git clone):
ENSEMBL from release 45

```
#From NCBI: 16 genomes
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_protein.faa.gz -o HOMSA_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_protein.faa.gz -o MUSMU_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_protein.faa.gz -o DANRE_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_protein.faa.gz -o GALGA_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.3_Xenopus_tropicalis_v9.1/GCF_000004195.3_Xenopus_tropicalis_v9.1_protein.faa.gz -o XENTR_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/165/045/GCF_000165045.1_Callorhinchus_milii-6.1.3/GCF_000165045.1_Callorhinchus_milii-6.1.3_protein.faa.gz -o CALMI_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/143/615/GCF_004143615.1_amil_sf_1.1/GCF_004143615.1_amil_sf_1.1_protein.faa.gz -o ACRMI_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/224/145/GCF_000224145.3_KH/GCF_000224145.3_KH_protein.faa.gz -o CIOIN_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/225/785/GCF_000225785.1_LatCha1/GCF_000225785.1_LatCha1_protein.faa.gz -o LATCH_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/955/GCF_000238955.4_M_zebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_protein.faa.gz -o MAYZE_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/242/695/GCF_000242695.1_LepOcu1/GCF_000242695.1_LepOcu1_protein.faa.gz -o LEPOC_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/344/595/GCF_000344595.1_CheMyd_1.0/GCF_000344595.1_CheMyd_1.0_protein.faa.gz -o CHEMY_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/365/465/GCF_000365465.2_Ptep_2.0/GCF_000365465.2_Ptep_2.0_protein.faa.gz -o PARTE_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/965/GCF_001417965.1_Aiptasia_genome_1.1/GCF_001417965.1_Aiptasia_genome_1.1_protein.faa.gz -o EXAPA_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/891/065/GCF_001891065.1_H_comes_QL1_v1/GCF_001891065.1_H_comes_QL1_v1_protein.faa.gz -o HIPCO_NCBI_raw_pep.fa.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/949/145/GCF_001949145.1_OKI-Apl_1.0/GCF_001949145.1_OKI-Apl_1.0_protein.faa.gz -o ACAPL_NCBI_raw_pep.fa.gz

#From ENSEMBL: 20 genomes
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/adineta_vaga/pep/Adineta_vaga.AMS_PRJEB1171_v1.pep.all.fa.gz -o ADIVA_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/amphimedon_queenslandica/pep/Amphimedon_queenslandica.Aqu1.pep.all.fa.gz -o AMPQU_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/anopheles_gambiae/pep/Anopheles_gambiae.AgamP4.pep.all.fa.gz -o ANOGA_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/release-45/metazoa/fasta/branchiostoma_lanceolatum/pep/Branchiostoma_lanceolatum.BraLan2.pep.all.fa.gz -o BRALA_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/caenorhabditis_elegans/pep/Caenorhabditis_elegans.WBcel235.pep.all.fa.gz -o CAEEL_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/capitella_teleta/pep/Capitella_teleta.Capitella_teleta_v1.0.pep.all.fa.gz -o CAPTE_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/crassostrea_gigas/pep/Crassostrea_gigas.oyster_v9.pep.all.fa.gz -o CRAGI_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/daphnia_pulex/pep/Daphnia_pulex.V1.0.pep.all.fa.gz -o DAPPU_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/drosophila_melanogaster/pep/Drosophila_melanogaster.BDGP6.22.pep.all.fa.gz -o DROME_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/helobdella_robusta/pep/Helobdella_robusta.Helro1.pep.all.fa.gz -o HELRO_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/ixodes_scapularis/pep/Ixodes_scapularis.IscaW1.pep.all.fa.gz -o IXOSC_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/lingula_anatina/pep/Lingula_anatina.LinAna1.0.pep.all.fa.gz -o LINAN_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/lottia_gigantea/pep/Lottia_gigantea.Lotgi1.pep.all.fa.gz -o LOTGI_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/nematostella_vectensis/pep/Nematostella_vectensis.ASM20922v1.pep.all.fa.gz -o NEMVE_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/strigamia_maritima/pep/Strigamia_maritima.Smar1.pep.all.fa.gz -o STRMA_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/strongylocentrotus_purpuratus/pep/Strongylocentrotus_purpuratus.Spur_3.1.pep.all.fa.gz -o STRPU_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/tribolium_castaneum/pep/Tribolium_castaneum.Tcas5.2.pep.all.fa.gz -o TRICA_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/trichoplax_adhaerens/pep/Trichoplax_adhaerens.ASM15027v1.pep.all.fa.gz -o TRIAD_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/protists/release-45/fasta/protists_choanoflagellida1_collection/salpingoeca_rosetta_gca_000188695/pep/Salpingoeca_rosetta_gca_000188695.Proterospongia_sp_ATCC50818.pep.all.fa.gz -o SALRO_ENSEMBL_raw_pep.fa.gz
curl ftp://ftp.ensemblgenomes.org/pub/protists/release-45/fasta/protists_ichthyosporea1_collection/capsaspora_owczarzaki_atcc_30864_gca_000151315/pep/Capsaspora_owczarzaki_atcc_30864_gca_000151315.C_owczarzaki_V2.pep.all.fa.gz -o CAPOW_ENSEMBL_raw_pep.fa.gz

#From other sources: 7 genomes
curl http://marimba.obs-vlfr.fr/download/file/fid/50 -o CLYHE_OTHER_raw_pep.fa
curl https://research.nhgri.nih.gov/hydra/download/genemodels_proteins/hydra2.0_genemodels.aa.gz -o HYDVU_OTHER_raw_pep.fa.gz
git clone https://bitbucket.org/molpalmuc/hoilungia-genome/src/master/sequences/Hhon_BRAKER1_proteins.fasta.gz HOIHO_raw_pep.fa.gz
cp /proj/Simakov/OTHER/SCALLOP/PYgenome-pep.fa MIZYE_OTHER_raw_pep.fa
curl https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.aa.gz -o MNELE_OTHER_raw_pep.fa.gz
curl https://marinegenomics.oist.jp/acornworm/download/pfl_public_ver1.0.prot > PTYFL_OTHER_raw_pep.fa
```
The 7 remaining genomes were downloaded with a browser (and or built using annotation files with gffread):

* Aurelia aurita from: https://drive.google.com/drive/folders/1NC6bZ9cxWkZyofOsMPzrxIH3C7m1ySiu, Aurelia.Genome_v1.2_Protein_Models_12-28-18.fasta, renamed AURAU_raw_pep.fa
* Euprymna scolopes clusters peptide file was copied from: /proj/Simakov/EUPRYMNA/MAHDITRAN/Euprymna_scolopes.fa , named EUPSC_raw_pep.fa
* Hofstenia miamia gff3 file was downloaded from http://srivastavalab.rc.fas.harvard.edu gffread was used to build HOFMI_raw_pep.fa using hmi_gene_annotation.gff3 and hmi_genome.fa
* Saccoglossus kowalevskii: SkowalevskiiJGIv3.0.longestTrs.pep.fa.gz Downloaded from Metazome v3 (https://metazome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Skowalevskii_er), called SACKO_OTHER_raw_pep.fa.gz
* Schmidtea mediterranea: gffread used to build SCHME_OTHER_ra_pep.fa using gff http://planmine.mpi-cbg.de/planmine/model/bulkdata/smes_v2_hconf_SMESG.gff3.zip (high confidence transcripts) and  http://planmine.mpi-cbg.de/planmine/model/bulkdata/dd_Smes_g4.fasta.zip


```
gffread SCHME_OTHER.gff3 -g dd_Smes_g4.fasta -y SCHME_OTHER_raw_pep.fa
```
* Sycon ciliatum genome, CDS and peptide downloaded from datadryad (https://datadryad.org/resource/doi:10.5061/dryad.tn0f3).
* Pleurobrachia bachei filtered gene models mRNAs (CDS) were downloaded from https://neurobase.rc.ufl.edu/pleurobrachia/download (03_P-bachei_Filtered_Gene_Models_RNA.txt) and also translated into peptides (also using transeq, frame 1), file was named PLEBA_raw_pep.fa
	we use
	```sed -i 's/_1$//' PLEBA.fa``` to delete the frame added by transseq to the accession

Finally, we gunzip all the gz files 

```
gunzip *.fa.gz
```



### 1.2. Download annotation files
### 1.2.1 Download of available gff files
gff files were available for 47 out of the 49 species. For SYCCI and PLEBA, we'll map the transcripts (see 1.2.2)
39/49 downloaded using curl:
We want the gffs which show all the sequences. For example, in the case of a chromosomal level assembly, we want the unplaced scaffolds as well.

```
#16 Gffs downloaded from NCBI
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/949/145/GCF_001949145.1_OKI-Apl_1.0/GCF_001949145.1_OKI-Apl_1.0_genomic.gff.gz -o ACAPL_NCBI.gff3.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/143/615/GCF_004143615.1_amil_sf_1.1/GCF_004143615.1_amil_sf_1.1_genomic.gff.gz -o ACRMI_NCBI.gff.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/165/045/GCF_000165045.1_Callorhinchus_milii-6.1.3/GCF_000165045.1_Callorhinchus_milii-6.1.3_genomic.gff.gz -o CALMI_NCBI.gff3.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/224/145/GCF_000224145.3_KH/GCF_000224145.3_KH_genomic.gff.gz -o CIOIN_NCBI.gff3.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.gff.gz -o DANRE_NCBI.gff3.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/965/GCF_001417965.1_Aiptasia_genome_1.1/GCF_001417965.1_Aiptasia_genome_1.1_genomic.gff.gz -o EXAPA_NCBI.gff.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.gff.gz -o GALGA_NCBI.gff3.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/242/695/GCF_000242695.1_LepOcu1/GCF_000242695.1_LepOcu1_genomic.gff.gz -o LEPOC_NCBI.gff3.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.gff.gz -o MUSMU_NCBI.gff3.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gff.gz -o HOMSA_NCBI.gff.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/365/465/GCF_000365465.2_Ptep_2.0/GCF_000365465.2_Ptep_2.0_genomic.gff.gz -o PARTE_NCBI.gff3.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/225/785/GCF_000225785.1_LatCha1/GCF_000225785.1_LatCha1_genomic.gff.gz -o LATCH_NCBI.gff.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/891/065/GCF_001891065.1_H_comes_QL1_v1/GCF_001891065.1_H_comes_QL1_v1_genomic.gff.gz -o HIPCO_NCBI.gff.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/955/GCF_000238955.4_M_zebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.gff.gz -o MAYZE_NCBI.gff.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/344/595/GCF_000344595.1_CheMyd_1.0/GCF_000344595.1_CheMyd_1.0_genomic.gff.gz -o CHEMY_NCBI.gff.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.3_Xenopus_tropicalis_v9.1/GCF_000004195.3_Xenopus_tropicalis_v9.1_genomic.gff.gz -o XENTR_NCBI.gff.gz

#20 Gffs downloaded from ENSEMBL
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/adineta_vaga/Adineta_vaga.AMS_PRJEB1171_v1.45.gff3.gz -o ADIVA_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/anopheles_gambiae//Anopheles_gambiae.AgamP4.45.gff3.gz -o ANOGA_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/caenorhabditis_elegans//Caenorhabditis_elegans.WBcel235.45.gff3.gz -o CAEEL_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/capitella_teleta//Capitella_teleta.Capitella_teleta_v1.0.45.gff3.gz -o CAPTE_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/crassostrea_gigas/Crassostrea_gigas.oyster_v9.45.gff3.gz -o CRAGI_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/daphnia_pulex/Daphnia_pulex.V1.0.45.gff3.gz -o DAPPU_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/drosophila_melanogaster/Drosophila_melanogaster.BDGP6.22.45.gff3.gz -o DROME_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/helobdella_robusta/Helobdella_robusta.Helro1.45.gff3.gz -o HELRO_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/ixodes_scapularis/Ixodes_scapularis.IscaW1.45.gff3.gz -o IXOSC_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/lingula_anatina/Lingula_anatina.LinAna1.0.45.gff3.gz -o LINAN_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/lottia_gigantea/Lottia_gigantea.Lotgi1.45.gff3.gz -o LOTGI_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/nematostella_vectensis//Nematostella_vectensis.ASM20922v1.45.gff3.gz -o NEMVE_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/strigamia_maritima//Strigamia_maritima.Smar1.45.gff3.gz -o STRMA_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/strongylocentrotus_purpuratus//Strongylocentrotus_purpuratus.Spur_3.1.45.gff3.gz -o STRPU_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/tribolium_castaneum/Tribolium_castaneum.Tcas5.2.45.gff3.gz -o TRICA_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/gff3/trichoplax_adhaerens//Trichoplax_adhaerens.ASM15027v1.45.gff3.gz -o TRIAD_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/protists/release-45/gff3/protists_choanoflagellida1_collection/salpingoeca_rosetta_gca_000188695/Salpingoeca_rosetta_gca_000188695.Proterospongia_sp_ATCC50818.45.gff3.gz -o SALRO_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/protists/release-45/gff3/protists_ichthyosporea1_collection/capsaspora_owczarzaki_atcc_30864_gca_000151315/Capsaspora_owczarzaki_atcc_30864_gca_000151315.C_owczarzaki_V2.45.gff3.gz -o CAPOW_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/release-45/metazoa/gff3/amphimedon_queenslandica/Amphimedon_queenslandica.Aqu1.45.gff3.gz -o AMPQU_ENSEMBL.gff3.gz
curl ftp://ftp.ensemblgenomes.org/pub/release-45/metazoa/gff3/branchiostoma_lanceolatum/Branchiostoma_lanceolatum.BraLan2.45.gff3.gz -o BRALA_ENSEMBL.gff3.gz

#4 gffs downloaded from other ressources
curl https://research.nhgri.nih.gov/hydra/download/genemodels_gff3/hydra2.0_genemodels.gff3.gz -o HYDVU_OTHER.gff3.gz
curl curl https://bitbucket.org/molpalmuc/hoilungia-genome/raw/0d523a5b8556741a37918f3f30d0ed0414833912/tracks/Hhon_BRAKER1_CDS.gff3.gz -o HOIHO_OTHER.gff3.gz
curl https://marinegenomics.oist.jp/acornworm/download/pfl_public_ver1.0.gff3 > PTYFL_OTHER.gff3
curl https://research.nhgri.nih.gov/mnemiopsis/download/proteome/ML2.2.gff3.gz -o MNELE_OTHER.gff3.gz
```

* *Aurelia aurita* (https://drive.google.com/drive/folders/1NC6bZ9cxWkZyofOsMPzrxIH3C7m1ySiu), named AURAU_OTHER.gff3
* *Branchiostoma lanceolatum*, gtf downloaded from https://www.dropbox.com/s/d4fqnoa8gdix3pa/Bla_annot_final.gtf.gz named BRALA_OTHER.gtf.gz
* *Clytia hemisphaerica*, gff from http://marimba.obs-vlfr.fr/download/file/fid/51 called CLYHE_OTHER.gff3
* *Euprymna scolopes*: gff located in /proj/Simakov/EUPRYMNA/MAHDITRAN/clusters_esc.gff3 use gffparser in the next parts firectly on this file.
* *Hofstenia miamia*, gff3 downloaded from http://srivastavalab.rc.fas.harvard.edu, hmi_gene_annotation.gff3 was renamed HOFMI_OTHER.gff3
* *Hoilungia hongkongensis*, gff3 downloaded from https://bitbucket.org/molpalmuc/hoilungia-genome/src/master/tracks/Hhon_BRAKER1_CDS.gff3.gz named HOIHO_OTHER.gff3
* *Mizuhopecten yessoensis* gff located in /proj/Simakov/OTHER/SCALLOP/chr.id.gff3 use gffparser in the next parts directly on this file.
* *Schmidtea mediterranea*: donwloaded high confidence gene predictions from  http://planmine.mpi-cbg.de/planmine/aspect.do?name=Gene%20Predictions, called the file SCHME_OTHER.gff3
* *Saccoglossus kowalevskii*: downloaded the longest transcrfipts form Metazome v3    SkowalevskiiJGIv3.0.longestTrs.gff3.gz

### 1.2.2 Mapping of transcripts/gene models onto the genome when annotation files were not uploaded by the authors (Pleurobrachia)
The annotations for *Pleurobrachia bachei* gene models and *Sycon ciliatum* transcripts were not available.
Genome database build using **gmap_gbuild**, transcripts mapped onto the genome using **gmap**

```
module load gmap
gmap_build -D Pleurobrachia_bachei -d PLEBA_genome GCA_000695325.1_P.bachei_draft_genome_v.1.1_genomic.fna
gmap -D Pleurobrachia_bachei/ -d PLEBA_genome 03_P-bachei_Filtered_Gene_Models_RNA.txt --gff3-add-separators=1 -f gff3_gene > PLEBA_03_mapped.gff3 

gmap_build -D Sycon_ciliatum/ -d SYCCI_genome sycon.genome.fa
gmap -D Sycon_ciliatum/ -d SYCCI_genome sycon.cds.fa --gff3-add-separators=1 -f gff3_gene > SYCCI_mapped.gff3 
```

* *Pleurobrachia bachei* : **18871 out of 18950** transcripts are mapped (i.e. 99.58%)
* *Sycon ciliatum* : **49278 out of 50731** transcripts are mapped (i.e. 97.00%)

We'll keep all the peptides in the following analyses (and not only the ones mapped onto the genome), because we want to differentiate between loss of genes and loss of synteny.
Keep the best mapped transcript
```
grep mrna1 SYCCI_mapped.gff3 | grep -P '\texon\t' > SYCCI_mapped_final.gff3
sed -i 's/sctid/scpid/g' SYCCI_mapped_final.gff3
grep mrna1 PLEBA_03_mapped.gff3 | grep -P '\texon\t' > PLEBA_mapped_final.gff3
```
## 2. Preparation of data for susbequent analyses
There is some preliminary steps we want to do before running the gff parser.

* For *H .miamia*, we modify the gff like this:

`sed -i '/\|/s/=.*\|/=/' HOFMI_hmi_gene_annotation.gff3`

* In the case of *S. ciliatum* and *P. bachei*, gmap outputs multiple paths per mRNA. We want to keep only one. Also, gmap maps both CDS and exons. Since we're mapping only CDS (some don't start with ATG), we'll keep only exons features.


```
grep mrna1 SYCCI_mapped.gff3 | grep -P '\texon\t' > SYCCI_mapped_final.gff3
sed -i 's/sctid/scpid/g' SYCCI_mapped_final.gff3
grep mrna1 PLEBA_03_mapped.gff3 | grep -P '\texon\t' > PLEBA_mapped_final.gff3
```

In the case of *S. kowalevskii* the accessions in the fasta files are gene names, and are not in the CDS lines. So we'll edit the gff with this script. It'll add up ";protein_id=Sakowv12345678m" comment in the last field of the gff on the lines containing CDS features. This is to make the file parseable by pyMakeMap

```
sk_gff_preparation.py SACKO_SkowalevskiiJGIv3.0.longestTrs.gff3 > SACKO_prepped.gff3
rm SACKO_SkowalevskiiJGIv3.0.longestTrs.gff3
```

### 2.1 Prepare non-redundant chrom and fasta files from gff
CHROM files are the format we commonly use. They're made by parsing gff annotation files. Some gffs comprise info about exon, gene, CDS, UTR.

We define gene coordinates as the CDS information for every species.

#### 2.1.1 Prepare chrom files for NCBI data
First, we'll delete the mitochondrial scaffolds from the gffs. Not all genome assemblies posess mitochondrial, we don't want mitochondrial cluster as falsely assigned synteny novelties.

```
grep "NC_007788.1" CALMI_NCBI.gff3| cut -f1 | grep -wvf - ACAPL_NCBI.gff3 > ACAPL_NCBI_filt.gff3; rm ACAPL_NCBI.gff3
#ACRMI has no MT chrom
grep "NC_014285.1" CALMI_NCBI.gff3| cut -f1 | grep -wvf - CALMI_NCBI.gff3 > CALMI_NCBI_filt.gff; rm CALMI_NCBI.gff3
grep "NC_017929.1" CIOIN_NCBI.gff3| cut -f1 | grep -wvf - CIOIN_NCBI.gff3 > CIOIN_NCBI_filt.gff3; rm CIOIN_NCBI.gff3
grep "NC_002333.2" DANRE_NCBI.gff3| cut -f1 | grep -wvf - DANRE_NCBI.gff3 > DANRE_filt_NCBI.gff3; rm DANRE_NCBI.gff3
#EXAPA has no MT chrom
grep -v "NC_040902.1" GALGA_NCBI.gff3 > GALGA_filt_NCBI.gff3; rm GALGA_NCBI.gff3
grep "NC_004744.1" LEPOC_NCBI.gff3| cut -f1 | grep -wvf - LEPOC_NCBI.gff3 > LEPOC_filt_NCBI.gff3; rm LEPOC_NCBI.gff3
grep "NC_005089.1" MUSMU_NCBI.gff3| cut -f1 | grep -wvf - MUSMU_NCBI.gff3 > MUSMU_filt_NCBI.gff3; rm MUSMU_NCBI.gff3
grep "NC_012920.1" HOMSA_NCBI.gff| cut -f1|   grep -wvf - HOMSA_NCBI.gff  > HOMSA_filt_NCBI.gff; rm HOMSA_NCBI.gff
#PARTE has no MT
grep "NC_001804.1" LATCH_NCBI.gff|cut -f1|grep -wvf - LATCH_NCBI.gff > LATCH_filt_NCBI.gff; rm LATCH_NCBI.gff
grep "NC_020336.1" HIPCO_NCBI.gff|cut -f1|grep -wvf - HIPCO_NCBI.gff > HIPCO_filt_NCBI.gff; rm HIPCO_NCBI.gff
grep "NC_027944.1" MAYZE_NCBI.gff|cut -f1|grep -wvf - MAYZE_NCBI.gff > MAYZE_filt_NCBI.gff; rm MAYZE_NCBI.gff
grep "NC_000886.1" CHEMY_NCBI.gff|cut -f1|grep -wvf - CHEMY_NCBI.gff > CHEMY_NCBI_filt.gff; rm CHEMY_NCBI.gff
grep "NC_006839.1" XENTR_NCBI.gff|cut -f1|grep -wvf - XENTR_NCBI.gff > XENTR_filt_NCBI.gff; rm XENTR_NCBI.gff
```

We use custom scripts to filters the longest isoforms (as identified by shared gene IDs). ENSEMBL gene ids are in the fasta headers. In the case of NCBI, gene IDs are not in the fasta file, so the annotation file (gff format) is also needed to determine to which gene the isoforms come from.

* ACAPL

```
LongestIsoforms_NCBI.py ../../proteins/ACAPL_NCBI_raw_pep.fa ACAPL_NCBI_filt.gff3 ACAPL_NCBI_longest_isoforms.fa
sed -i 's/^>/>ACAPL_/' ACAPL_NCBI_longest_isoforms.fa
pymakeMap.py -gff ACAPL_NCBI_filt.gff3 -p ACAPL -f CDS -k protein_id -o ../chrom/ACAPL.chrom -r fasta -F ACAPL_NCBI_longest_isoforms.fa
cut -f2 ../chrom/ACAPL.chrom > ACAPL.list ; pyfasta extract --file ACAPL.list --header --space --fasta ACAPL_NCBI_longest_isoforms.fa > ../../proteins_processed/ACAPL.fasta
rm ACAPL_NCBI_longest_isoforms* ACAPL.list
```

* ACRMI

```
LongestIsoforms_NCBI.py ../../proteins/ACRMI_NCBI_raw_pep.fa ACRMI_NCBI.gff ACRMI_NCBI_longest_isoforms.fa
sed -i 's/^>/>ACRMI_/' ACRMI_NCBI_longest_isoforms.fa
pymakeMap.py -gff ACRMI_NCBI.gff -p ACRMI -f CDS -k protein_id -o ../chrom/ACRMI.chrom -r fasta -F ACRMI_NCBI_longest_isoforms.fa
cut -f2 ../chrom/ACRMI.chrom > ACRMI.list ; pyfasta extract --file ACRMI.list --header --space --fasta ACRMI_NCBI_longest_isoforms.fa > ../../proteins_processed/ACRMI.fasta
rm ACRMI_NCBI_longest_isoforms* ACRMI.list
```

* CALMI

```
LongestIsoforms_NCBI.py ../../proteins/CALMI_NCBI_raw_pep.fa CALMI_NCBI.gff3 CALMI_NCBI_longest_isoforms.fa
sed -i 's/^>/>CALMI_/' CALMI_NCBI_longest_isoforms.fa
pymakeMap.py -gff CALMI_NCBI.gff3 -p CALMI -f CDS -k protein_id -o ../chrom/CALMI.chrom -r fasta -F CALMI_NCBI_longest_isoforms.fa
cut -f2 ../chrom/CALMI.chrom > CALMI.list ; pyfasta extract --file CALMI.list --header --space --fasta CALMI_NCBI_longest_isoforms.fa > ../../proteins_processed/CALMI.fasta
rm CALMI_NCBI_longest_isoforms* CALMI.list
```

* CHEMY

```
LongestIsoforms_NCBI.py ../../proteins/CHEMY_NCBI_raw_pep.fa CHEMY_NCBI_filt.gff CHEMY_NCBI_longest_isoforms.fa
sed -i 's/^>/>CHEMY_/' CHEMY_NCBI_longest_isoforms.fa
pymakeMap.py -gff CHEMY_NCBI_filt.gff -p CHEMY -f CDS -k protein_id -o ../chrom/CHEMY.chrom -r fasta -F CHEMY_NCBI_longest_isoforms.fa
cut -f2 ../chrom/CHEMY.chrom > CHEMY.list ; pyfasta extract --file CHEMY.list --header --space --fasta CHEMY_NCBI_longest_isoforms.fa > ../../proteins_processed/CHEMY.fasta
rm CHEMY_NCBI_longest_isoforms* CHEMY.list
```

* CIOIN

```
LongestIsoforms_NCBI.py ../../proteins/CIOIN_NCBI_raw_pep.fa CIOIN_NCBI_filt.gff3 CIOIN_NCBI_longest_isoforms.fa
sed -i 's/^>/>CIOIN_/' CIOIN_NCBI_longest_isoforms.fa
pymakeMap.py -gff CIOIN_NCBI_filt.gff3 -p CIOIN -f CDS -k protein_id -o ../chrom/CIOIN.chrom -r fasta -F CIOIN_NCBI_longest_isoforms.fa
cut -f2 ../chrom/CIOIN.chrom > CIOIN.list ; pyfasta extract --file CIOIN.list --header --space --fasta CIOIN_NCBI_longest_isoforms.fa > ../../proteins_processed/CIOIN.fasta
rm CIOIN_NCBI_longest_isoforms* CIOIN.list
```

* DANRE

```
LongestIsoforms_NCBI.py ../../proteins/DANRE_NCBI_raw_pep.fa DANRE_filt_NCBI.gff3 DANRE_NCBI_longest_isoforms.fa
sed -i 's/^>/>DANRE_/' DANRE_NCBI_longest_isoforms.fa
pymakeMap.py -gff DANRE_filt_NCBI.gff3 -p DANRE -f CDS -k protein_id -o ../chrom/DANRE.chrom -r fasta -F DANRE_NCBI_longest_isoforms.fa
cut -f2 ../chrom/DANRE.chrom > DANRE.list ; pyfasta extract --file DANRE.list --header --space --fasta DANRE_NCBI_longest_isoforms.fa > ../../proteins_processed/DANRE.fasta
rm DANRE_NCBI_longest_isoforms* DANRE.list
```

* EXAPA

```
LongestIsoforms_NCBI.py ../../proteins/EXAPA_NCBI_raw_pep.fa EXAPA_NCBI.gff EXAPA_NCBI_longest_isoforms.fa
sed -i 's/^>/>EXAPA_/' EXAPA_NCBI_longest_isoforms.fa
pymakeMap.py -gff EXAPA_NCBI.gff -p EXAPA -f CDS -k protein_id -o ../chrom/EXAPA.chrom -r fasta -F EXAPA_NCBI_longest_isoforms.fa
cut -f2 ../chrom/EXAPA.chrom > EXAPA.list ; pyfasta extract --file EXAPA.list --header --space --fasta EXAPA_NCBI_longest_isoforms.fa > ../../proteins_processed/EXAPA.fasta
rm EXAPA_NCBI_longest_isoforms* EXAPA.list
```

* GALGA

```
LongestIsoforms_NCBI.py ../../proteins/GALGA_NCBI_raw_pep.fa GALGA_filt_NCBI.gff3 GALGA_NCBI_longest_isoforms.fa
sed -i 's/^>/>GALGA_/' GALGA_NCBI_longest_isoforms.fa
pymakeMap.py -gff GALGA_filt_NCBI.gff3 -p GALGA -f CDS -k protein_id -o ../chrom/GALGA.chrom -r fasta -F GALGA_NCBI_longest_isoforms.fa
cut -f2 ../chrom/GALGA.chrom > GALGA.list ; pyfasta extract --file GALGA.list --header --space --fasta GALGA_NCBI_longest_isoforms.fa > ../../proteins_processed/GALGA.fasta
rm GALGA_NCBI_longest_isoforms* GALGA.list
```

* HIPCO

```
LongestIsoforms_NCBI.py ../../proteins/HIPCO_NCBI_raw_pep.fa HIPCO_filt_NCBI.gff HIPCO_NCBI_longest_isoforms.fa
sed -i 's/^>/>HIPCO_/' HIPCO_NCBI_longest_isoforms.fa
pymakeMap.py -gff HIPCO_filt_NCBI.gff -p HIPCO -f CDS -k protein_id -o ../chrom/HIPCO.chrom -r fasta -F HIPCO_NCBI_longest_isoforms.fa
cut -f2 ../chrom/HIPCO.chrom > HIPCO.list ; pyfasta extract --file HIPCO.list --header --space --fasta HIPCO_NCBI_longest_isoforms.fa > ../../proteins_processed/HIPCO.fasta
rm HIPCO_NCBI_longest_isoforms* HIPCO.list
```

* HOMSA

```
LongestIsoforms_NCBI.py ../../proteins/HOMSA_NCBI_raw_pep.fa HOMSA_filt_NCBI.gff HOMSA_NCBI_longest_isoforms.fa
sed -i 's/^>/>HOMSA_/' HOMSA_NCBI_longest_isoforms.fa
pymakeMap.py -gff HOMSA_filt_NCBI.gff -p HOMSA -f CDS -k protein_id -o ../chrom/HOMSA.chrom -r fasta -F HOMSA_NCBI_longest_isoforms.fa
cut -f2 ../chrom/HOMSA.chrom > HOMSA.list ; pyfasta extract --file HOMSA.list --header --space --fasta HOMSA_NCBI_longest_isoforms.fa > ../../proteins_processed/HOMSA.fasta
rm HOMSA_NCBI_longest_isoforms* HOMSA.list
```

* LATCH

```
LongestIsoforms_NCBI.py ../../proteins/LATCH_NCBI_raw_pep.fa LATCH_filt_NCBI.gff LATCH_NCBI_longest_isoforms.fa
sed -i 's/^>/>LATCH_/' LATCH_NCBI_longest_isoforms.fa
pymakeMap.py -gff LATCH_filt_NCBI.gff -p LATCH -f CDS -k protein_id -o ../chrom/LATCH.chrom -r fasta -F LATCH_NCBI_longest_isoforms.fa
cut -f2 ../chrom/LATCH.chrom > LATCH.list ; pyfasta extract --file LATCH.list --header --space --fasta LATCH_NCBI_longest_isoforms.fa > ../../proteins_processed/LATCH.fasta
rm LATCH_NCBI_longest_isoforms* LATCH.list
```

* LEPOC

```
LongestIsoforms_NCBI.py ../../proteins/LEPOC_NCBI_raw_pep.fa LEPOC_filt_NCBI.gff3 LEPOC_NCBI_longest_isoforms.fa
sed -i 's/^>/>LEPOC_/' LEPOC_NCBI_longest_isoforms.fa
pymakeMap.py -gff LEPOC_filt_NCBI.gff3 -p LEPOC -f CDS -k protein_id -o ../chrom/LEPOC.chrom -r fasta -F LEPOC_NCBI_longest_isoforms.fa
cut -f2 ../chrom/LEPOC.chrom > LEPOC.list ; pyfasta extract --file LEPOC.list --header --space --fasta LEPOC_NCBI_longest_isoforms.fa > ../../proteins_processed/LEPOC.fasta
rm LEPOC_NCBI_longest_isoforms* LEPOC.list
```

* MAYZE

```
LongestIsoforms_NCBI.py ../../proteins/MAYZE_NCBI_raw_pep.fa MAYZE_filt_NCBI.gff MAYZE_NCBI_longest_isoforms.fa
sed -i 's/^>/>MAYZE_/' MAYZE_NCBI_longest_isoforms.fa
pymakeMap.py -gff MAYZE_filt_NCBI.gff -p MAYZE -f CDS -k protein_id -o ../chrom/MAYZE.chrom -r fasta -F MAYZE_NCBI_longest_isoforms.fa
cut -f2 ../chrom/MAYZE.chrom > MAYZE.list ; pyfasta extract --file MAYZE.list --header --space --fasta MAYZE_NCBI_longest_isoforms.fa > ../../proteins_processed/MAYZE.fasta
rm MAYZE_NCBI_longest_isoforms* MAYZE.list
```

* MUSMU

```
LongestIsoforms_NCBI.py ../../proteins/MUSMU_NCBI_raw_pep.fa MUSMU_filt_NCBI.gff3 MUSMU_NCBI_longest_isoforms.fa
sed -i 's/^>/>MUSMU_/' MUSMU_NCBI_longest_isoforms.fa
pymakeMap.py -gff MUSMU_filt_NCBI.gff3 -p MUSMU -f CDS -k protein_id -o ../chrom/MUSMU.chrom -r fasta -F MUSMU_NCBI_longest_isoforms.fa
cut -f2 ../chrom/MUSMU.chrom > MUSMU.list ; pyfasta extract --file MUSMU.list --header --space --fasta MUSMU_NCBI_longest_isoforms.fa > ../../proteins_processed/MUSMU.fasta
rm MUSMU_NCBI_longest_isoforms* MUSMU.list
```

* PARTE

```
LongestIsoforms_NCBI.py ../../proteins/PARTE_NCBI_raw_pep.fa PARTE_NCBI.gff3 PARTE_NCBI_longest_isoforms.fa
sed -i 's/^>/>PARTE_/' PARTE_NCBI_longest_isoforms.fa
pymakeMap.py -gff PARTE_NCBI.gff3 -p PARTE -f CDS -k protein_id -o ../chrom/PARTE.chrom -r fasta -F PARTE_NCBI_longest_isoforms.fa
cut -f2 ../chrom/PARTE.chrom > PARTE.list ; pyfasta extract --file PARTE.list --header --space --fasta PARTE_NCBI_longest_isoforms.fa > ../../proteins_processed/PARTE.fasta
rm PARTE_NCBI_longest_isoforms* PARTE.list
```

* XENTR

```
LongestIsoforms_NCBI.py ../../proteins/XENTR_NCBI_raw_pep.fa XENTR_filt_NCBI.gff XENTR_NCBI_longest_isoforms.fa
sed -i 's/^>/>XENTR_/' XENTR_NCBI_longest_isoforms.fa
pymakeMap.py -gff XENTR_filt_NCBI.gff -p XENTR -f CDS -k protein_id -o ../chrom/XENTR.chrom -r fasta -F XENTR_NCBI_longest_isoforms.fa
cut -f2 ../chrom/XENTR.chrom > XENTR.list ; pyfasta extract --file XENTR.list --header --space --fasta XENTR_NCBI_longest_isoforms.fa > ../../proteins_processed/XENTR.fasta
rm XENTR_NCBI_longest_isoforms* XENTR.list
```

#### 2.1.2 Prepare chrom files for ENSEMBL data

* ADIVA

```
LongestIsoforms_ENSEMBL.py ../../proteins/ADIVA_ENSEMBL_raw_pep.fa ADIVA_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>ADIVA_/' ADIVA_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff ADIVA_ENSEMBL.gff3 -p ADIVA -f CDS -k protein_id -o ../chrom/ADIVA.chrom -r fasta -F ADIVA_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/ADIVA.chrom > ADIVA.list ; pyfasta extract --file ADIVA.list --header --space --fasta ADIVA_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/ADIVA.fasta
rm ADIVA_ENSEMBL_longest_isoforms* ADIVA.list
```

* AMPQU

```
LongestIsoforms_ENSEMBL.py ../../proteins/AMPQU_ENSEMBL_raw_pep.fa AMPQU_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>AMPQU_/' AMPQU_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff AMPQU_ENSEMBL.gff3 -p AMPQU -f CDS -k protein_id -o ../chrom/AMPQU.chrom -r fasta -F AMPQU_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/AMPQU.chrom > AMPQU.list ; pyfasta extract --file AMPQU.list --header --space --fasta AMPQU_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/AMPQU.fasta
rm AMPQU_ENSEMBL_longest_isoforms* AMPQU.list
```

* ANOGA

```
grep -vP "Mt\tVectorBase"  ANOGA_ENSEMBL.gff3 > ANOGA_filt_ENSEMBL.gff3; rm ANOGA_ENSEMBL.gff3
LongestIsoforms_ENSEMBL.py ../../proteins/ANOGA_ENSEMBL_raw_pep.fa ANOGA_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>ANOGA_/' ANOGA_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff ANOGA_filt_ENSEMBL.gff3 -p ANOGA -f CDS -k protein_id -o ../chrom/ANOGA.chrom -r fasta -F ANOGA_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/ANOGA.chrom > ANOGA.list ; pyfasta extract --file ANOGA.list --header --space --fasta ANOGA_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/ANOGA.fasta
rm ANOGA_ENSEMBL_longest_isoforms* ANOGA.list
```

* BRALA

```
LongestIsoforms_ENSEMBL.py ../../proteins/BRALA_ENSEMBL_raw_pep.fa BRALA_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>BRALA_/' BRALA_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff BRALA_ENSEMBL.gff3 -p BRALA -f CDS -k protein_id -o ../chrom/BRALA.chrom -r fasta -F BRALA_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/BRALA.chrom > BRALA.list ; pyfasta extract --file BRALA.list --header --space --fasta BRALA_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/BRALA.fasta
rm BRALA_ENSEMBL_longest_isoforms* BRALA.list
```

* CAEEL

```
grep -vP "MtDNA\tWormBase"  CAEEL_ENSEMBL.gff3 > CAEEL_filt_ENSEMBL.gff3; rm CAEEL_ENSEMBL.gff3
LongestIsoforms_ENSEMBL.py ../../proteins/CAEEL_ENSEMBL_raw_pep.fa CAEEL_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>CAEEL_/' CAEEL_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff CAEEL_filt_ENSEMBL.gff3 -p CAEEL -f CDS -k protein_id -o ../chrom/CAEEL.chrom -r fasta -F CAEEL_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/CAEEL.chrom > CAEEL.list ; pyfasta extract --file CAEEL.list --header --space --fasta CAEEL_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/CAEEL.fasta
rm CAEEL_ENSEMBL_longest_isoforms* CAEEL.list
```

* CAPOW

```
LongestIsoforms_ENSEMBL.py ../../proteins/CAPOW_ENSEMBL_raw_pep.fa CAPOW_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>CAPOW_/' CAPOW_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff CAPOW_ENSEMBL.gff3 -p CAPOW -f CDS -k protein_id -o ../chrom/CAPOW.chrom -r fasta -F CAPOW_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/CAPOW.chrom > CAPOW.list ; pyfasta extract --file CAPOW.list --header --space --fasta CAPOW_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/CAPOW.fasta
rm CAPOW_ENSEMBL_longest_isoforms* CAPOW.list
```

* CAPTE

```
LongestIsoforms_ENSEMBL.py ../../proteins/CAPTE_ENSEMBL_raw_pep.fa CAPTE_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>CAPTE_/' CAPTE_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff CAPTE_ENSEMBL.gff3 -p CAPTE -f CDS -k protein_id -o ../chrom/CAPTE.chrom -r fasta -F CAPTE_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/CAPTE.chrom > CAPTE.list ; pyfasta extract --file CAPTE.list --header --space --fasta CAPTE_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/CAPTE.fasta
rm CAPTE_ENSEMBL_longest_isoforms* CAPTE.list
```

* CRAGI

```
LongestIsoforms_ENSEMBL.py ../../proteins/CRAGI_ENSEMBL_raw_pep.fa CRAGI_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>CRAGI_/' CRAGI_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff CRAGI_ENSEMBL.gff3 -p CRAGI -f CDS -k protein_id -o ../chrom/CRAGI.chrom -r fasta -F CRAGI_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/CRAGI.chrom > CRAGI.list ; pyfasta extract --file CRAGI.list --header --space --fasta CRAGI_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/CRAGI.fasta
rm CRAGI_ENSEMBL_longest_isoforms* CRAGI.list
```

* DAPPU

```
LongestIsoforms_ENSEMBL.py ../../proteins/DAPPU_ENSEMBL_raw_pep.fa DAPPU_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>DAPPU_/' DAPPU_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff DAPPU_ENSEMBL.gff3 -p DAPPU -f CDS -k protein_id -o ../chrom/DAPPU.chrom -r fasta -F DAPPU_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/DAPPU.chrom > DAPPU.list ; pyfasta extract --file DAPPU.list --header --space --fasta DAPPU_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/DAPPU.fasta
rm DAPPU_ENSEMBL_longest_isoforms* DAPPU.list
```

* DROME

```
grep -vP "mitochondrion_genome\tFlyBase"  DROME_ENSEMBL.gff3 > DROME_filt_ENSEMBL.gff3; rm DROME_ENSEMBL.gff3
LongestIsoforms_ENSEMBL.py ../../proteins/DROME_ENSEMBL_raw_pep.fa DROME_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>DROME_/' DROME_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff DROME_filt_ENSEMBL.gff3 -p DROME -f CDS -k protein_id -o ../chrom/DROME.chrom -r fasta -F DROME_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/DROME.chrom > DROME.list ; pyfasta extract --file DROME.list --header --space --fasta DROME_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/DROME.fasta
rm DROME_ENSEMBL_longest_isoforms* DROME.list
```

* HELRO

```
LongestIsoforms_ENSEMBL.py ../../proteins/HELRO_ENSEMBL_raw_pep.fa HELRO_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>HELRO_/' HELRO_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff HELRO_ENSEMBL.gff3 -p HELRO -f CDS -k protein_id -o ../chrom/HELRO.chrom -r fasta -F HELRO_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/HELRO.chrom > HELRO.list ; pyfasta extract --file HELRO.list --header --space --fasta HELRO_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/HELRO.fasta
rm HELRO_ENSEMBL_longest_isoforms* HELRO.list
```

* IXOSC

```
LongestIsoforms_ENSEMBL.py ../../proteins/IXOSC_ENSEMBL_raw_pep.fa IXOSC_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>IXOSC_/' IXOSC_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff IXOSC_ENSEMBL.gff3 -p IXOSC -f CDS -k protein_id -o ../chrom/IXOSC.chrom -r fasta -F IXOSC_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/IXOSC.chrom > IXOSC.list ; pyfasta extract --file IXOSC.list --header --space --fasta IXOSC_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/IXOSC.fasta
rm IXOSC_ENSEMBL_longest_isoforms* IXOSC.list
```

* LINAN

```
LongestIsoforms_ENSEMBL.py ../../proteins/LINAN_ENSEMBL_raw_pep.fa LINAN_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>LINAN_/' LINAN_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff LINAN_ENSEMBL.gff3 -p LINAN -f CDS -k protein_id -o ../chrom/LINAN.chrom -r fasta -F LINAN_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/LINAN.chrom > LINAN.list ; pyfasta extract --file LINAN.list --header --space --fasta LINAN_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/LINAN.fasta
rm LINAN_ENSEMBL_longest_isoforms* LINAN.list
```

* LOTGI

```
LongestIsoforms_ENSEMBL.py ../../proteins/LOTGI_ENSEMBL_raw_pep.fa LOTGI_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>LOTGI_/' LOTGI_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff LOTGI_ENSEMBL.gff3 -p LOTGI -f CDS -k protein_id -o ../chrom/LOTGI.chrom -r fasta -F LOTGI_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/LOTGI.chrom > LOTGI.list ; pyfasta extract --file LOTGI.list --header --space --fasta LOTGI_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/LOTGI.fasta
rm LOTGI_ENSEMBL_longest_isoforms* LOTGI.list
```

* NEMVE

```
LongestIsoforms_ENSEMBL.py ../../proteins/NEMVE_ENSEMBL_raw_pep.fa NEMVE_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>NEMVE_/' NEMVE_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff NEMVE_ENSEMBL.gff3 -p NEMVE -f CDS -k protein_id -o ../chrom/NEMVE.chrom -r fasta -F NEMVE_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/NEMVE.chrom > NEMVE.list ; pyfasta extract --file NEMVE.list --header --space --fasta NEMVE_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/NEMVE.fasta
rm NEMVE_ENSEMBL_longest_isoforms* NEMVE.list
```

* SALRO

```
LongestIsoforms_ENSEMBL.py ../../proteins/SALRO_ENSEMBL_raw_pep.fa SALRO_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>SALRO_/' SALRO_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff SALRO_ENSEMBL.gff3 -p SALRO -f CDS -k protein_id -o ../chrom/SALRO.chrom -r fasta -F SALRO_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/SALRO.chrom > SALRO.list ; pyfasta extract --file SALRO.list --header --space --fasta SALRO_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/SALRO.fasta
rm SALRO_ENSEMBL_longest_isoforms* SALRO.list
```

* STRMA

```
LongestIsoforms_ENSEMBL.py ../../proteins/STRMA_ENSEMBL_raw_pep.fa STRMA_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>STRMA_/' STRMA_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff STRMA_ENSEMBL.gff3 -p STRMA -f CDS -k protein_id -o ../chrom/STRMA.chrom -r fasta -F STRMA_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/STRMA.chrom > STRMA.list ; pyfasta extract --file STRMA.list --header --space --fasta STRMA_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/STRMA.fasta
rm STRMA_ENSEMBL_longest_isoforms* STRMA.list
```

* STRPU, we'll use `-d` flag of pymakemap to delete isoforms registered with different gene ids. It deletes 2891 proteins
This flag deletes overlapping genes with the following criteria:
* isoform with 1 exon is merged with overlapping isoform of 1+ size if any exon is overlapping on more than 20 base pairs
* isoform with 2 exons is merged with overlapping isoform of 2+ size if exons share at least a start or a stop
* isoform with 3 exons is merged with overlapping isoform of 3+ size if they share at least one internal exon


```
LongestIsoforms_ENSEMBL.py ../../proteins/STRPU_ENSEMBL_raw_pep.fa STRPU_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>STRPU_/' STRPU_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff STRPU_ENSEMBL.gff3 -p STRPU -f CDS -k protein_id -o ../chrom/STRPU.chrom -r fasta -F STRPU_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/STRPU.chrom > STRPU.list ; pyfasta extract --file STRPU.list --header --space --fasta STRPU_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/STRPU.fasta
rm STRPU_ENSEMBL_longest_isoforms* STRPU.list
```

* TRIAD

```
LongestIsoforms_ENSEMBL.py ../../proteins/TRIAD_ENSEMBL_raw_pep.fa TRIAD_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>TRIAD_/' TRIAD_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff TRIAD_ENSEMBL.gff3 -p TRIAD -f CDS -k protein_id -o ../chrom/TRIAD.chrom -r fasta -F TRIAD_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/TRIAD.chrom > TRIAD.list ; pyfasta extract --file TRIAD.list --header --space --fasta TRIAD_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/TRIAD.fasta
rm TRIAD_ENSEMBL_longest_isoforms* TRIAD.list
```

* TRICA

```
LongestIsoforms_ENSEMBL.py ../../proteins/TRICA_ENSEMBL_raw_pep.fa TRICA_ENSEMBL_longest_isoforms.fa
sed -i 's/^>/>TRICA_/' TRICA_ENSEMBL_longest_isoforms.fa
pymakeMap.py -gff TRICA_ENSEMBL.gff3 -p TRICA -f CDS -k protein_id -o ../chrom/TRICA.chrom -r fasta -F TRICA_ENSEMBL_longest_isoforms.fa
cut -f2 ../chrom/TRICA.chrom > TRICA.list ; pyfasta extract --file TRICA.list --header --space --fasta TRICA_ENSEMBL_longest_isoforms.fa > ../../proteins_processed/TRICA.fasta
rm TRICA_ENSEMBL_longest_isoforms* TRICA.list
```

#### 2.1.3 Prepare chrom files for the data from OTHER sources
##### 2.1.3.1 Prepare data with transcript variants
* HOIHO

```
LongestIsoforms_HOIHO_HYDVU_PTYFL.py ../../proteins/HOIHO_OTHER_raw_pep.fa HOIHO_OTHER_longest_isoforms.fa
sed -i 's/^>/>HOIHO_/' HOIHO_OTHER_longest_isoforms.fa
pymakeMap.py -gff HOIHO_OTHER.gff3 -p HOIHO -f CDS -k Name -o ../chrom/HOIHO.chrom -r fasta -F HOIHO_OTHER_longest_isoforms.fa
cut -f2 ../chrom/HOIHO.chrom > HOIHO.list ; pyfasta extract --file HOIHO.list --header --space --fasta HOIHO_OTHER_longest_isoforms.fa > ../../proteins_processed/HOIHO.fasta
rm HOIHO_OTHER_longest_isoforms* HOIHO.list
```

* HYDVU

```
LongestIsoforms_HOIHO_HYDVU_PTYFL.py ../../proteins/HYDVU_OTHER_raw_pep.fa HYDVU_OTHER_longest_isoforms.fa
sed -i 's/^>/>HYDVU_/' HYDVU_OTHER_longest_isoforms.fa
pymakeMap.py -gff HYDVU_OTHER.gff3 -p HYDVU -f CDS -k Parent -o ../chrom/HYDVU.chrom -r fasta -F HYDVU_OTHER_longest_isoforms.fa
cut -f2 ../chrom/HYDVU.chrom > HYDVU.list ; pyfasta extract --file HYDVU.list --header --space --fasta HYDVU_OTHER_longest_isoforms.fa > ../../proteins_processed/HYDVU.fasta
rm HYDVU_OTHER_longest_isoforms* HYDVU.list
```

* PTYFL

```
LongestIsoforms_HOIHO_HYDVU_PTYFL.py ../../proteins/PTYFL_OTHER_raw_pep.fa PTYFL_OTHER_longest_isoforms.fa
sed -i 's/^>/>PTYFL_/' PTYFL_OTHER_longest_isoforms.fa
pymakeMap.py -gff PTYFL_OTHER.gff3 -p PTYFL -f CDS -k Parent -o ../chrom/PTYFL.chrom -r fasta -F PTYFL_OTHER_longest_isoforms.fa
cut -f2 ../chrom/PTYFL.chrom > PTYFL.list ; pyfasta extract --file PTYFL.list --header --space --fasta PTYFL_OTHER_longest_isoforms.fa > ../../proteins_processed/PTYFL.fasta
rm PTYFL_OTHER_longest_isoforms* PTYFL.list
```

* SCHME

```
LongestIsoforms_SCHME.py ../../proteins/SCHME_OTHER_raw_pep.fa SCHME_OTHER.gff3 SCHME_OTHER_longest_isoforms.fa
sed -i 's/^>/>SCHME_/' SCHME_OTHER_longest_isoforms.fa
pymakeMap.py -gff SCHME_OTHER.gff3 -p SCHME -f CDS -k Parent -o ../chrom/SCHME.chrom -r fasta -F SCHME_OTHER_longest_isoforms.fa
cut -f2 ../chrom/SCHME.chrom > SCHME.list ; pyfasta extract --file SCHME.list --header --space --fasta SCHME_OTHER_longest_isoforms.fa > ../../proteins_processed/SCHME.fasta
rm SCHME_OTHER_longest_isoforms* SCHME.list
```

##### 2.1.3.2 Prepare data without transcript variants
For all the following animals, proteomes are already filtered for redundancy

* AURAU

```
sed 's/^>/>AURAU_/' ../../proteins/AURAU_OTHER_raw_pep.fa > AURAU_OTHER_longest_isoforms.fa
pymakeMap.py -gff AURAU_other.gff3 -p AURAU -f CDS -k protein_id -o ../chrom/AURAU.chrom -r fasta -F AURAU_OTHER_longest_isoforms.fa
cut -f2 ../chrom/AURAU.chrom > AURAU.list ; pyfasta extract --file AURAU.list --header --space --fasta AURAU_OTHER_longest_isoforms.fa > ../../proteins_processed/AURAU.fasta
rm AURAU_OTHER_longest_isoforms* AURAU.list
```

* CLYHE

```
sed 's/^>\([^-]*\)-protein/>CLYHE_\1/' ../../proteins/CLYHE_OTHER_raw_pep.fa > CLYHE_OTHER_longest_isoforms.fa
pymakeMap.py -gff CLYHE_OTHER.gff3 -p CLYHE -f exon -k Parent -o ../chrom/CLYHE.chrom -r fasta -F CLYHE_OTHER_longest_isoforms.fa
cut -f2 ../chrom/CLYHE.chrom > CLYHE.list ; pyfasta extract --file CLYHE.list --header --space --fasta CLYHE_OTHER_longest_isoforms.fa > ../../proteins_processed/CLYHE.fasta
rm CLYHE_OTHER_longest_isoforms* CLYHE.list
```

* EUPSC

```
sed 's/\tcluster_/\tID=cluster_/' esc_allchroms_uniq_15jul2019.gff > EUPSC.gff
pymakeMap.py -gff EUPSC.gff -p EUPSC -f exon -k ID -o ../chrom/EUPSC.chrom -r fasta -F EUPSC_OTHER_longest_isoforms.fa
cut -f2 ../chrom/EUPSC.chrom > EUPSC.list ; pyfasta extract --file EUPSC.list --header --space --fasta EUPSC_OTHER_longest_isoforms.fa > ../../proteins_processed/EUPSC.fasta
rm EUPSC_OTHER_longest_isoforms* EUPSC.list
```

* HOFMI

```
sed 's/^>\([^_]\+_[^_]\+\)_\(g[0-9]\+.t[0-9]\)/>HOFMI_\2\ \1/' ../../proteins/HOFMI_OTHER_raw_pep.fa > HOFMI_OTHER_longest_isoforms.fa
pymakeMap.py -gff HOFMI_OTHER.gff3 -p HOFMI -f CDS -k Parent -o ../chrom/HOFMI.chrom -r fasta -F HOFMI_OTHER_longest_isoforms.fa
cut -f2 ../chrom/HOFMI.chrom > HOFMI.list ; pyfasta extract --file HOFMI.list --header --space --fasta HOFMI_OTHER_longest_isoforms.fa > ../../proteins_processed/HOFMI.fasta
rm HOFMI_OTHER_longest_isoforms* HOFMI.list
```

* MIZYE

```
sed 's/^>/>MIZYE_/' ../../proteins/MIZYE_OTHER_raw_pep.fa > MIZYE_OTHER_longest_isoforms.fa
pymakeMap.py -gff SCALLOP/chr.id.gff3 -p MIZYE -f CDS -k Parent -o ../chrom/MIZYE.chrom -r fasta -F MIZYE_OTHER_longest_isoforms.fa
cut -f2 ../chrom/MIZYE.chrom > MIZYE.list ; pyfasta extract --file MIZYE.list --header --space --fasta MIZYE_OTHER_longest_isoforms.fa > ../../proteins_processed/MIZYE.fasta
rm MIZYE_OTHER_longest_isoforms* MIZYE.list
```

* MNELE

```
sed 's/^>/>MNELE_/' ../../proteins/MNELE_OTHER_raw_pep.fa > MNELE_OTHER_longest_isoforms.fa
pymakeMap.py -gff MNELE_OTHER.gff3 -p MNELE -f CDS -k Parent -o ../chrom/MNELE.chrom -r fasta -F MNELE_OTHER_longest_isoforms.fa
cut -f2 ../chrom/MNELE.chrom > MNELE.list ; pyfasta extract --file MNELE.list --header --space --fasta MNELE_OTHER_longest_isoforms.fa > ../../proteins_processed/MNELE.fasta
rm MNELE_OTHER_longest_isoforms* MNELE.list
```

* PLEBA

```
sed 's/^>/>PLEBA_/' ../../proteins/PLEBA_OTHER_raw_pep.fa > PLEBA_OTHER_longest_isoforms.fa
pymakeMap.py -gff PLEBA_mapped_final.gff3 -p PLEBA -f exon -k Name -o ../chrom/PLEBA.chrom -r fasta -F PLEBA_OTHER_longest_isoforms.fa
cut -f2 ../chrom/PLEBA.chrom > PLEBA.list ; pyfasta extract --file PLEBA.list --header --space --fasta PLEBA_OTHER_longest_isoforms.fa > ../../proteins_processed/PLEBA.fasta
rm PLEBA_OTHER_longest_isoforms* PLEBA.list
```

* SACKO

```
sed 's/^>/>SACKO_/' ../../proteins/SACKO_OTHER_raw_pep.fa > SACKO_OTHER_longest_isoforms.fa
pymakeMap.py -gff SACKO_prepped.gff3 -p SACKO -f CDS -k protein_id -o ../chrom/SACKO.chrom -r fasta -F SACKO_OTHER_longest_isoforms.fa
cut -f2 ../chrom/SACKO.chrom > SACKO.list ; pyfasta extract --file SACKO.list --header --space --fasta SACKO_OTHER_longest_isoforms.fa > ../../proteins_processed/SACKO.fasta
rm SACKO_OTHER_longest_isoforms* SACKO.list
```

* SYCCI

```
sed 's/^>/>SYCCI_/' ../../proteins/SYCCI_OTHER_raw_pep.fa > SYCCI_OTHER_longest_isoforms.fa
pymakeMap.py -gff SYCCI_mapped_final.gff3 -p SYCCI -f exon -k Name -o ../chrom/SYCCI.chrom -r fasta -F SYCCI_OTHER_longest_isoforms.fa -d
cut -f2 ../chrom/SYCCI.chrom > SYCCI.list ; pyfasta extract --file SYCCI.list --header --space --fasta SYCCI_OTHER_longest_isoforms.fa > ../../proteins_processed/SYCCI.fasta
rm SYCCI_OTHER_longest_isoforms* SYCCI.list
```
## 3. Orthology assignment

First, we want to generate the blast commands in preparation for orthofinder.
Small check that the proteins do not end with symbols such as * or . to represent translated STOP codons, this trips up makeblastdb.

```
sed -i 's/*//' EUPSC.fasta
sed -i 's/*//' PLEBA.fasta
sed -i 's/*//' SACKO.fasta
sed -i 's/\([A-Z]\)\.\([A-Z]\)/\1\2/g' HOFMI.fasta
sed -i 's/\([A-Z]\)\.\([A-Z]\)/\1\2/g' SCHME.fasta

module load orthofinder ncbiblastplus
orthofinder -S blast -f 00_orthology_assignment/proteins_processed/ -op
```


Next, we want to enqueue the jobs. We thus use slurmtasks and tmprewrite to build the slurm script (https://github.com/nijibabulu/slurm-utils).

```
module load ncbiblastplus


for i in $(seq 0 48); do touch 00_orthology_assignment/proteins_processed/OrthoFinder/Results_Nov15/WorkingDirectory/BlastDBSpecies$i; done
for i in $(seq 0 48); do
  for j in $(seq 0 48); do
    faspeciesi="00_orthology_assignment/proteins_processed/OrthoFinder/Results_Nov15/WorkingDirectory/Species${i}.fa"
    DBspeciesj="00_orthology_assignment/proteins_processed/OrthoFinder/Results_Nov15/WorkingDirectory/BlastDBSpecies${j}"
    outfile="00_orthology_assignment/proteins_processed/OrthoFinder/Results_Nov15/WorkingDirectory/Blast${i}_${j}.txt"
    tmprewrite "echo {$DBspeciesj.phr:i} {$DBspeciesj.pin:i} {$DBspeciesj.psq:i}; blastp -num_threads 8 -outfmt 6 -evalue 0.001 -query {$faspeciesi:i} -db {$DBspeciesj:i} -out {$outfile:o}"
  done
done > 00_orthology_assignment/proteins_processed/ofblast
sed -i 's/echo\ [^;]*;\ //' ofblast

slurmtasks -p 8core 00_orthology_assignment/proteins_processed/ofblast| sbatch --array=1-192
slurmtasks -p 8core 00_orthology_assignment/proteins_processed/ofblast| sbatch
```

Now we want to run the graph building (Orthofinder) and clustering (MCL).
```
module load orthofinder
orthofinder -og -b 00_orthology_assignment/proteins_processed/OrthoFinder/Results_Nov15/WorkingDirectory
# convert file to clus
orthoFinderToOrthogroup.pl 00_orthology_assignment/proteins_processed/OrthoFinder/Results_Nov15/WorkingDirectory/Orthogroups.txt > 01_microsynteny/Orthofinder.clus
```## 4. Infer microsyntenic blocks
The two perl scripts `prepMicroSynt.pl` and `makeClusters3.pl` can be found at https://github.com/nijibabulu/metazoan_synteny/tree/master/scripts. These scripts are part of the microsynteny pipeline intially published in Simakov et al. 2013  (10.1038/nature11696):


### 4.1 randomize genomes
We'll randomize the genomes 3 times. Using the randomizations for both approaches.
we use `ls` to get a list of chrom files, and pymakeRandChrom takes the list from stdin.
```
ls 01_microsynteny/chrom/*.chrom | pymakeRandChrom.py -n 3 -o 01_microsynteny/randomized_chrom
```

### 4.2 Microsynteny computed using Orthofinder's OGs
First, let's compute the microsynteny with the Orthofinder clus file

```
cd 01_microsynteny/chrom
prepMicroSynt.pl ACAPL.chrom,ACRMI.chrom,ADIVA.chrom,AMPQU.chrom,ANOGA.chrom,AURAU.chrom,BRALA.chrom,CAEEL.chrom,CALMI.chrom,CAPOW.chrom,CAPTE.chrom,CHEMY.chrom,CIOIN.chrom,CLYHE.chrom,CRAGI.chrom,DANRE.chrom,DAPPU.chrom,DROME.chrom,EUPSC.chrom,EXAPA.chrom,GALGA.chrom,HELRO.chrom,HIPCO.chrom,HOFMI.chrom,HOIHO.chrom,HOMSA.chrom,HYDVU.chrom,IXOSC.chrom,LATCH.chrom,LEPOC.chrom,LINAN.chrom,LOTGI.chrom,MAYZE.chrom,MIZYE.chrom,MNELE.chrom,MUSMU.chrom,NEMVE.chrom,PARTE.chrom,PLEBA.chrom,PTYFL.chrom,SACKO.chrom,SALRO.chrom,SCHME.chrom,STRMA.chrom,STRPU.chrom,SYCCI.chrom,TRIAD.chrom,TRICA.chrom,XENTR.chrom 5 01_microsynteny/Orthofinder.clus
sbatch --array=1-1176 --constraint=array-1core --job-name=synt_of_OBS job.sh

#reorganize a bit the files.
mkdir bash_scripts log_files pairwise_blocks; mv *.out log_files; mv *.sh bash_scripts; mv *.blocks pairwise_blocks

CHROMS=../../chrom/ACAPL.chrom,../../chrom/ACRMI.chrom,../../chrom/ADIVA.chrom,../../chrom/AMPQU.chrom,../../chrom/ANOGA.chrom,../../chrom/AURAU.chrom,../../chrom/BRALA.chrom,../../chrom/CAEEL.chrom,../../chrom/CALMI.chrom,../../chrom/CAPOW.chrom,../../chrom/CAPTE.chrom,../../chrom/CHEMY.chrom,../../chrom/CIOIN.chrom,../../chrom/CLYHE.chrom,../../chrom/CRAGI.chrom,../../chrom/DANRE.chrom,../../chrom/DAPPU.chrom,../../chrom/DROME.chrom,../../chrom/EUPSC.chrom,../../chrom/EXAPA.chrom,../../chrom/GALGA.chrom,../../chrom/HELRO.chrom,../../chrom/HIPCO.chrom,../../chrom/HOFMI.chrom,../../chrom/HOIHO.chrom,../../chrom/HOMSA.chrom,../../chrom/HYDVU.chrom,../../chrom/IXOSC.chrom,../../chrom/LATCH.chrom,../../chrom/LEPOC.chrom,../../chrom/LINAN.chrom,../../chrom/LOTGI.chrom,../../chrom/MAYZE.chrom,../../chrom/MIZYE.chrom,../../chrom/MNELE.chrom,../../chrom/MUSMU.chrom,../../chrom/NEMVE.chrom,../../chrom/PARTE.chrom,../../chrom/PLEBA.chrom,../../chrom/PTYFL.chrom,../../chrom/SACKO.chrom,../../chrom/SALRO.chrom,../../chrom/SCHME.chrom,../../chrom/STRMA.chrom,../../chrom/STRPU.chrom,../../chrom/SYCCI.chrom,../../chrom/TRIAD.chrom,../../chrom/TRICA.chrom,../../chrom/XENTR.chrom

cd 01_microsynteny/chrom_of/pairwise_blocks
makeClusters3.pl $CHROMS .5.blocks 3 0.3 0.5 &> 5.blocks.3.syn.synt

#correct the coordinates in the synt file, coordinates outputted by perl scripts are not exact block boundaries.
correct_blocks_coordinates.py 5.blocks.3.syn.synt $CHROMS > 5.blocks.3.syn_corrected.synt
```

Now, for the first randomized genome

```
cd 01_microsynteny/randomized_chrom/rand.1/Orthofinder
prepMicroSynt.pl ../ACAPL.chrom.rand.1,../ACRMI.chrom.rand.1,../ADIVA.chrom.rand.1,../AMPQU.chrom.rand.1,../ANOGA.chrom.rand.1,../AURAU.chrom.rand.1,../BRALA.chrom.rand.1,../CAEEL.chrom.rand.1,../CALMI.chrom.rand.1,../CAPOW.chrom.rand.1,../CAPTE.chrom.rand.1,../CHEMY.chrom.rand.1,../CIOIN.chrom.rand.1,../CLYHE.chrom.rand.1,../CRAGI.chrom.rand.1,../DANRE.chrom.rand.1,../DAPPU.chrom.rand.1,../DROME.chrom.rand.1,../EUPSC.chrom.rand.1,../EXAPA.chrom.rand.1,../GALGA.chrom.rand.1,../HELRO.chrom.rand.1,../HIPCO.chrom.rand.1,../HOFMI.chrom.rand.1,../HOIHO.chrom.rand.1,../HOMSA.chrom.rand.1,../HYDVU.chrom.rand.1,../IXOSC.chrom.rand.1,../LATCH.chrom.rand.1,../LEPOC.chrom.rand.1,../LINAN.chrom.rand.1,../LOTGI.chrom.rand.1,../MAYZE.chrom.rand.1,../MIZYE.chrom.rand.1,../MNELE.chrom.rand.1,../MUSMU.chrom.rand.1,../NEMVE.chrom.rand.1,../PARTE.chrom.rand.1,../PLEBA.chrom.rand.1,../PTYFL.chrom.rand.1,../SACKO.chrom.rand.1,../SALRO.chrom.rand.1,../SCHME.chrom.rand.1,../STRMA.chrom.rand.1,../STRPU.chrom.rand.1,../SYCCI.chrom.rand.1,../TRIAD.chrom.rand.1,../TRICA.chrom.rand.1,../XENTR.chrom.rand.1 5 01_microsynteny/Orthofinder.clus
sbatch --array=1-1176 --constraint=array-1core --job-name=synt_of_RAND1 job.sh

CHROMS_RAND1="../../ACAPL.chrom.rand.1,../../ACRMI.chrom.rand.1,../../ADIVA.chrom.rand.1,../../AMPQU.chrom.rand.1,../../ANOGA.chrom.rand.1,../../AURAU.chrom.rand.1,../../BRALA.chrom.rand.1,../../CAEEL.chrom.rand.1,../../CALMI.chrom.rand.1,../../CAPOW.chrom.rand.1,../../CAPTE.chrom.rand.1,../../CHEMY.chrom.rand.1,../../CIOIN.chrom.rand.1,../../CLYHE.chrom.rand.1,../../CRAGI.chrom.rand.1,../../DANRE.chrom.rand.1,../../DAPPU.chrom.rand.1,../../DROME.chrom.rand.1,../../EUPSC.chrom.rand.1,../../EXAPA.chrom.rand.1,../../GALGA.chrom.rand.1,../../HELRO.chrom.rand.1,../../HIPCO.chrom.rand.1,../../HOFMI.chrom.rand.1,../../HOIHO.chrom.rand.1,../../HOMSA.chrom.rand.1,../../HYDVU.chrom.rand.1,../../IXOSC.chrom.rand.1,../../LATCH.chrom.rand.1,../../LEPOC.chrom.rand.1,../../LINAN.chrom.rand.1,../../LOTGI.chrom.rand.1,../../MAYZE.chrom.rand.1,../../MIZYE.chrom.rand.1,../../MNELE.chrom.rand.1,../../MUSMU.chrom.rand.1,../../NEMVE.chrom.rand.1,../../PARTE.chrom.rand.1,../../PLEBA.chrom.rand.1,../../PTYFL.chrom.rand.1,../../SACKO.chrom.rand.1,../../SALRO.chrom.rand.1,../../SCHME.chrom.rand.1,../../STRMA.chrom.rand.1,../../STRPU.chrom.rand.1,../../SYCCI.chrom.rand.1,../../TRIAD.chrom.rand.1,../../TRICA.chrom.rand.1,../../XENTR.chrom.rand.1"

mkdir bash_scripts log_files pairwise_blocks; mv *.out log_files; mv *.sh bash_scripts; mv *.blocks pairwise_blocks
cd 01_microsynteny/randomized_chrom/rand.1/Orthofinder/pairwise_blocks
makeClusters3.pl $CHROMS_RAND1 .5.blocks 3 0.3 0.5 > 5.blocks.3.syn.synt

correct_blocks_coordinates.py 5.blocks.3.syn.synt $CHROMS_RAND1 > 5.blocks.3.syn_corrected.synt
```

The second randomized genome

```
cd 01_microsynteny/randomized_chrom/rand.2/Orthofinder
prepMicroSynt.pl ../ACAPL.chrom.rand.2,../ACRMI.chrom.rand.2,../ADIVA.chrom.rand.2,../AMPQU.chrom.rand.2,../ANOGA.chrom.rand.2,../AURAU.chrom.rand.2,../BRALA.chrom.rand.2,../CAEEL.chrom.rand.2,../CALMI.chrom.rand.2,../CAPOW.chrom.rand.2,../CAPTE.chrom.rand.2,../CHEMY.chrom.rand.2,../CIOIN.chrom.rand.2,../CLYHE.chrom.rand.2,../CRAGI.chrom.rand.2,../DANRE.chrom.rand.2,../DAPPU.chrom.rand.2,../DROME.chrom.rand.2,../EUPSC.chrom.rand.2,../EXAPA.chrom.rand.2,../GALGA.chrom.rand.2,../HELRO.chrom.rand.2,../HIPCO.chrom.rand.2,../HOFMI.chrom.rand.2,../HOIHO.chrom.rand.2,../HOMSA.chrom.rand.2,../HYDVU.chrom.rand.2,../IXOSC.chrom.rand.2,../LATCH.chrom.rand.2,../LEPOC.chrom.rand.2,../LINAN.chrom.rand.2,../LOTGI.chrom.rand.2,../MAYZE.chrom.rand.2,../MIZYE.chrom.rand.2,../MNELE.chrom.rand.2,../MUSMU.chrom.rand.2,../NEMVE.chrom.rand.2,../PARTE.chrom.rand.2,../PLEBA.chrom.rand.2,../PTYFL.chrom.rand.2,../SACKO.chrom.rand.2,../SALRO.chrom.rand.2,../SCHME.chrom.rand.2,../STRMA.chrom.rand.2,../STRPU.chrom.rand.2,../SYCCI.chrom.rand.2,../TRIAD.chrom.rand.2,../TRICA.chrom.rand.2,../XENTR.chrom.rand.2 5 01_microsynteny/Orthofinder.clus
sbatch --array=1-1176 --constraint=array-1core --job-name=synt_of_RAND2 job.sh

CHROM_RAND2="../../ACAPL.chrom.rand.2,../../ACRMI.chrom.rand.2,../../ADIVA.chrom.rand.2,../../AMPQU.chrom.rand.2,../../ANOGA.chrom.rand.2,../../AURAU.chrom.rand.2,../../BRALA.chrom.rand.2,../../CAEEL.chrom.rand.2,../../CALMI.chrom.rand.2,../../CAPOW.chrom.rand.2,../../CAPTE.chrom.rand.2,../../CHEMY.chrom.rand.2,../../CIOIN.chrom.rand.2,../../CLYHE.chrom.rand.2,../../CRAGI.chrom.rand.2,../../DANRE.chrom.rand.2,../../DAPPU.chrom.rand.2,../../DROME.chrom.rand.2,../../EUPSC.chrom.rand.2,../../EXAPA.chrom.rand.2,../../GALGA.chrom.rand.2,../../HELRO.chrom.rand.2,../../HIPCO.chrom.rand.2,../../HOFMI.chrom.rand.2,../../HOIHO.chrom.rand.2,../../HOMSA.chrom.rand.2,../../HYDVU.chrom.rand.2,../../IXOSC.chrom.rand.2,../../LATCH.chrom.rand.2,../../LEPOC.chrom.rand.2,../../LINAN.chrom.rand.2,../../LOTGI.chrom.rand.2,../../MAYZE.chrom.rand.2,../../MIZYE.chrom.rand.2,../../MNELE.chrom.rand.2,../../MUSMU.chrom.rand.2,../../NEMVE.chrom.rand.2,../../PARTE.chrom.rand.2,../../PLEBA.chrom.rand.2,../../PTYFL.chrom.rand.2,../../SACKO.chrom.rand.2,../../SALRO.chrom.rand.2,../../SCHME.chrom.rand.2,../../STRMA.chrom.rand.2,../../STRPU.chrom.rand.2,../../SYCCI.chrom.rand.2,../../TRIAD.chrom.rand.2,../../TRICA.chrom.rand.2,../../XENTR.chrom.rand.2"

mkdir bash_scripts log_files pairwise_blocks; mv *.out log_files; mv *.sh bash_scripts; mv *.blocks pairwise_blocks
cd 01_microsynteny/randomized_chrom/rand.1/Orthofinder/pairwise_blocks
makeClusters3.pl $CHROM_RAND2 3 0.3 0.5 > 5.blocks.3.syn.synt

correct_blocks_coordinates.py 5.blocks.3.syn.synt $CHROM_RAND2 > 5.blocks.3.syn_corrected.synt
```

And the third

```
cd 01_microsynteny/randomized_chrom/rand.3/Orthofinder
prepMicroSynt.pl ../ACAPL.chrom.rand.3,../ACRMI.chrom.rand.3,../ADIVA.chrom.rand.3,../AMPQU.chrom.rand.3,../ANOGA.chrom.rand.3,../AURAU.chrom.rand.3,../BRALA.chrom.rand.3,../CAEEL.chrom.rand.3,../CALMI.chrom.rand.3,../CAPOW.chrom.rand.3,../CAPTE.chrom.rand.3,../CHEMY.chrom.rand.3,../CIOIN.chrom.rand.3,../CLYHE.chrom.rand.3,../CRAGI.chrom.rand.3,../DANRE.chrom.rand.3,../DAPPU.chrom.rand.3,../DROME.chrom.rand.3,../EUPSC.chrom.rand.3,../EXAPA.chrom.rand.3,../GALGA.chrom.rand.3,../HELRO.chrom.rand.3,../HIPCO.chrom.rand.3,../HOFMI.chrom.rand.3,../HOIHO.chrom.rand.3,../HOMSA.chrom.rand.3,../HYDVU.chrom.rand.3,../IXOSC.chrom.rand.3,../LATCH.chrom.rand.3,../LEPOC.chrom.rand.3,../LINAN.chrom.rand.3,../LOTGI.chrom.rand.3,../MAYZE.chrom.rand.3,../MIZYE.chrom.rand.3,../MNELE.chrom.rand.3,../MUSMU.chrom.rand.3,../NEMVE.chrom.rand.3,../PARTE.chrom.rand.3,../PLEBA.chrom.rand.3,../PTYFL.chrom.rand.3,../SACKO.chrom.rand.3,../SALRO.chrom.rand.3,../SCHME.chrom.rand.3,../STRMA.chrom.rand.3,../STRPU.chrom.rand.3,../SYCCI.chrom.rand.3,../TRIAD.chrom.rand.3,../TRICA.chrom.rand.3,../XENTR.chrom.rand.3 5 01_microsynteny/Orthofinder.clus
sbatch --array=1-1176 --constraint=array-1core --job-name=synt_of_RAND3 job.sh

CHROM_RAND3=../../ACAPL.chrom.rand.3,../../ACRMI.chrom.rand.3,../../ADIVA.chrom.rand.3,../../AMPQU.chrom.rand.3,../../ANOGA.chrom.rand.3,../../AURAU.chrom.rand.3,../../BRALA.chrom.rand.3,../../CAEEL.chrom.rand.3,../../CALMI.chrom.rand.3,../../CAPOW.chrom.rand.3,../../CAPTE.chrom.rand.3,../../CHEMY.chrom.rand.3,../../CIOIN.chrom.rand.3,../../CLYHE.chrom.rand.3,../../CRAGI.chrom.rand.3,../../DANRE.chrom.rand.3,../../DAPPU.chrom.rand.3,../../DROME.chrom.rand.3,../../EUPSC.chrom.rand.3,../../EXAPA.chrom.rand.3,../../GALGA.chrom.rand.3,../../HELRO.chrom.rand.3,../../HIPCO.chrom.rand.3,../../HOFMI.chrom.rand.3,../../HOIHO.chrom.rand.3,../../HOMSA.chrom.rand.3,../../HYDVU.chrom.rand.3,../../IXOSC.chrom.rand.3,../../LATCH.chrom.rand.3,../../LEPOC.chrom.rand.3,../../LINAN.chrom.rand.3,../../LOTGI.chrom.rand.3,../../MAYZE.chrom.rand.3,../../MIZYE.chrom.rand.3,../../MNELE.chrom.rand.3,../../MUSMU.chrom.rand.3,../../NEMVE.chrom.rand.3,../../PARTE.chrom.rand.3,../../PLEBA.chrom.rand.3,../../PTYFL.chrom.rand.3,../../SACKO.chrom.rand.3,../../SALRO.chrom.rand.3,../../SCHME.chrom.rand.3,../../STRMA.chrom.rand.3,../../STRPU.chrom.rand.3,../../SYCCI.chrom.rand.3,../../TRIAD.chrom.rand.3,../../TRICA.chrom.rand.3,../../XENTR.chrom.rand.3
mkdir bash_scripts log_files pairwise_blocks; mv *.out log_files; mv *.sh bash_scripts; mv *.blocks pairwise_blocks

cd 01_microsynteny/randomized_chrom/rand.3/Orthofinder/pairwise_blocks
makeClusters3.pl $CHROM_RAND3 .5.blocks 3 0.3 0.5 > 5.blocks.3.syn.synt

correct_blocks_coordinates.py 5.blocks.3.syn.synt $CHROM_RAND3 > 5.blocks.3.syn_corrected.synt
```
## 5. Gene density analysis
### 5.1 Data preparation
#### 5.1.1 Download genome files
If possible, sequences have been downloaded as softmasked (repetitive sequences in lowercase). This is more informative than hardmask (all residues in uppercase, masked residues replaces by Ns).

```
#From NCBI: 16 genomes
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.39_GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.fna.gz -o HOMSA.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/635/GCF_000001635.26_GRCm38.p6/GCF_000001635.26_GRCm38.p6_genomic.fna.gz -o MUSMU.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.fna.gz -o DANRE.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_genomic.fna.gz -o GALGA.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.3_Xenopus_tropicalis_v9.1/GCF_000004195.3_Xenopus_tropicalis_v9.1_genomic.fna.gz -o XENTR.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/165/045/GCF_000165045.1_Callorhinchus_milii-6.1.3/GCF_000165045.1_Callorhinchus_milii-6.1.3_genomic.fna.gz -o CALMI.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/004/143/615/GCF_004143615.1_amil_sf_1.1/GCF_004143615.1_amil_sf_1.1_genomic.fna.gz -o ACRMI.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/224/145/GCF_000224145.3_KH/GCF_000224145.3_KH_genomic.fna.gz -o CIOIN.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/225/785/GCF_000225785.1_LatCha1/GCF_000225785.1_LatCha1_genomic.fna.gz -o LATCH.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/238/955/GCF_000238955.4_M_zebra_UMD2a/GCF_000238955.4_M_zebra_UMD2a_genomic.fna.gz -o MAYZE.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/242/695/GCF_000242695.1_LepOcu1/GCF_000242695.1_LepOcu1_genomic.fna.gz -o LEPOC.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/344/595/GCF_000344595.1_CheMyd_1.0/GCF_000344595.1_CheMyd_1.0_genomic.fna.gz -o CHEMY.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/365/465/GCF_000365465.2_Ptep_2.0/GCF_000365465.2_Ptep_2.0_genomic.fna.gz -o PARTE.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/417/965/GCF_001417965.1_Aiptasia_genome_1.1/GCF_001417965.1_Aiptasia_genome_1.1_genomic.fna.gz -o EXAPA.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/891/065/GCF_001891065.1_H_comes_QL1_v1/GCF_001891065.1_H_comes_QL1_v1_genomic.fna.gz -o HIPCO.masked_genome.gz
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/949/145/GCF_001949145.1_OKI-Apl_1.0/GCF_001949145.1_OKI-Apl_1.0_genomic.fna.gz -o ACAPL.masked_genome.gz

#From ENSEMBL: 20 genomes
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/adineta_vaga/dna/Adineta_vaga.AMS_PRJEB1171_v1.dna_sm.toplevel.fa.gz -o ADIVA.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/amphimedon_queenslandica/dna/Amphimedon_queenslandica.Aqu1.dna_sm.toplevel.fa.gz -o AMPQU.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/anopheles_gambiae/dna/Anopheles_gambiae.AgamP4.dna_sm.toplevel.fa.gz -o ANOGA.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/release-45/metazoa/fasta/branchiostoma_lanceolatum/dna/Branchiostoma_lanceolatum.BraLan2.dna_sm.toplevel.fa.gz -o BRALA.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/caenorhabditis_elegans/dna/Caenorhabditis_elegans.WBcel235.dna_sm.toplevel.fa.gz -o CAEEL.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/capitella_teleta/dna/Capitella_teleta.Capitella_teleta_v1.0.dna_sm.toplevel.fa.gz -o CAPTE.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/crassostrea_gigas/dna/Crassostrea_gigas.oyster_v9.dna_sm.toplevel.fa.gz -o CRAGI.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/daphnia_pulex/dna/Daphnia_pulex.V1.0.dna_sm.toplevel.fa.gz -o DAPPU.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP6.22.dna_sm.toplevel.fa.gz -o DROME.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/helobdella_robusta/dna/Helobdella_robusta.Helro1.dna_sm.toplevel.fa.gz -o HELRO.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/ixodes_scapularis/dna/Ixodes_scapularis.IscaW1.dna_sm.toplevel.fa.gz -o IXOSC.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/lingula_anatina/dna/Lingula_anatina.LinAna1.0.dna_sm.toplevel.fa.gz -o LINAN.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/lottia_gigantea/dna/Lottia_gigantea.Lotgi1.dna_sm.toplevel.fa.gz -o LOTGI.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/nematostella_vectensis/dna/Nematostella_vectensis.ASM20922v1.dna_sm.toplevel.fa.gz -o NEMVE.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/strigamia_maritima/dna/Strigamia_maritima.Smar1.dna_sm.toplevel.fa.gz -o STRMA.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/strongylocentrotus_purpuratus/dna/Strongylocentrotus_purpuratus.Spur_3.1.dna_sm.toplevel.fa.gz -o STRPU.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/tribolium_castaneum/dna/Tribolium_castaneum.Tcas5.2.dna_sm.toplevel.fa.gz -o TRICA.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/metazoa/release-45/fasta/trichoplax_adhaerens/dna/Trichoplax_adhaerens.ASM15027v1.dna_sm.toplevel.fa.gz -o TRIAD.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/protists/release-45/fasta/protists_choanoflagellida1_collection/salpingoeca_rosetta_gca_000188695/dna/Salpingoeca_rosetta_gca_000188695.Proterospongia_sp_ATCC50818.dna_sm.toplevel.fa.gz -o SALRO.masked_genome.gz
curl ftp://ftp.ensemblgenomes.org/pub/protists/release-45/fasta/protists_ichthyosporea1_collection/capsaspora_owczarzaki_atcc_30864_gca_000151315/dna/Capsaspora_owczarzaki_atcc_30864_gca_000151315.C_owczarzaki_V2.dna_sm.toplevel.fa.gz -o CAPOW.masked_genome.gz

#From other sources: 7 genomes
curl http://marimba.obs-vlfr.fr/download/file/fid/48 -o CLYHE.hardmasked_genome.gz
curl https://research.nhgri.nih.gov/hydra/download/assembly/Hm105_Dovetail_Assembly_1.0.fa.gz -o HYDVU.masked_genome.gz
curl https://bitbucket.org/molpalmuc/hoilungia-genome/raw/0d523a5b8556741a37918f3f30d0ed0414833912/sequences/Hhon_final_contigs_softmasked.fasta.gz -o HOIHO.masked_genome.gz
curl https://research.nhgri.nih.gov/mnemiopsis/download/genome/MlScaffold09.nt.gz -o MNELE.masked_genome.gz
curl https://marinegenomics.oist.jp/acornworm/download/pfl_scaffold_ver1.0.14.masked.fasta.gz -o PTYFL.masked_genome.gz
curl http://planmine.mpi-cbg.de/planmine/model/bulkdata/dd_Smes_g4.fasta.zip -o SCHME.masked_genome
curl ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/695/325/GCA_000695325.1_P.bachei_draft_genome_v.1.1/GCA_000695325.1_P.bachei_draft_genome_v.1.1_genomic.fna.gz -o PLEBA.masked_genome
```

The 3 remaining genomes were downloaded with a browser (and or made using annotation files) or used from another emplacment in cube:

* Aurelia aurita from: https://drive.google.com/drive/folders/1NC6bZ9cxWkZyofOsMPzrxIH3C7m1ySiu, Aurelia.Genome_v1.2_Protein_Models_12-28-18.fasta, renamed AURAU.hardmasked_genome
* Saccoglossus kowalevskii: SkowalevskiiJGIv3.0.longestTrs.pep.fa.gz Downloaded from Metazome v3 (https://metazome.jgi.doe.gov/pz/portal.html#!bulk?org=Org_Skowalevskii_er), called SACKO.masked_genome.gz
* Sycon ciliatum genome, CDS and peptide downloaded from datadryad (https://datadryad.org/resource/doi:10.5061/dryad.tn0f3).


EUPSC.masked_genome genome from Schmidbaur et al. in prep.

Finally, we gunzip all the gz files 

```
gunzip *.fa.gz
```
## 5 Gene density analysis (REDUX)
### 5.1. reconstructing ancestral counts of syntenic blocks
We'll use 5 alternating topologies:

* Ctenophore-sister (CS) (cteno,(sponge,(plac,(cnid,(acoel,prot,deut)))));
* Sponge-sister (PS) (sponge,(cteno,(plac,(cnid,(acoel,prot,deut)))));
* Xenacoelomorpha (XNEP): (cteno,sponge,(plac,(cnid,(acoel(prot,deut)))))); 
* Xenambulacraria (XAMB): (cteno,sponge,(plac,(cnid,(prot((chord,(acoel,ambulacraria)))))));
* Consensus (CONS) (cteno,sponge,(plac,(cnid,(acoel,prot,deut))));

Here are the actual trees:

* CONS
```
(((((((((((((((HOMSA,MUSMU),(CHEMY,GALGA)),XENTR),LATCH),(((MAYZE,HIPCO),DANRE),LEPOC)),CALMI)Vertebrata,CIOIN)Olfactores,BRALA)Chordata,((SACKO,PTYFL),(STRPU,ACAPL))Ambulacraria)Deuterostomia,(((EUPSC,LOTGI,(MIZYE,CRAGI))Mollusca,(CAPTE,HELRO),ADIVA,LINAN,SCHME)Lophotrochozoa,((((((DROME,ANOGA),TRICA),DAPPU),STRMA),(IXOSC,PARTE)),CAEEL)Ecdysozoa)Protostomia,HOFMI)Bilateria,(((NEMVE,EXAPA),ACRMI),((HYDVU,CLYHE),AURAU))Cnidaria)Planulozoa,(HOIHO,TRIAD)Placozoa)Parahoxozoa,(PLEBA,MNELE)Ctenophora,(SYCCI,AMPQU)Porifera)Metazoa,SALRO)Choanozoa,CAPOW)Filozoa;
```
* CS
```
((((((((((((((((HOMSA,MUSMU),(CHEMY,GALGA)),XENTR),LATCH),(((MAYZE,HIPCO),DANRE),LEPOC)),CALMI)Vertebrata,CIOIN)Olfactores,BRALA)Chordata,((SACKO,PTYFL),(STRPU,ACAPL))Ambulacraria)Deuterostomia,(((EUPSC,LOTGI,(MIZYE,CRAGI))Mollusca,(CAPTE,HELRO),ADIVA,LINAN,SCHME)Lophotrochozoa,((((((DROME,ANOGA),TRICA),DAPPU),STRMA),(IXOSC,PARTE)),CAEEL)Ecdysozoa)Protostomia,HOFMI)Bilateria,(((NEMVE,EXAPA),ACRMI),((HYDVU,CLYHE),AURAU))Cnidaria)Planulozoa,(HOIHO,TRIAD)Placozoa)Parahoxozoa,(SYCCI,AMPQU)Porifera)Porifera_Parahoxozoa,(PLEBA,MNELE)Ctenophora)Metazoa,SALRO)Choanozoa,CAPOW)Filozoa;
```
* PS
```
((((((((((((((((HOMSA,MUSMU),(CHEMY,GALGA)),XENTR),LATCH),(((MAYZE,HIPCO),DANRE),LEPOC)),CALMI)Vertebrata,CIOIN)Olfactores,BRALA)Chordata,((SACKO,PTYFL),(STRPU,ACAPL))Ambulacraria)Deuterostomia,(((EUPSC,LOTGI,(MIZYE,CRAGI))Mollusca,(CAPTE,HELRO),ADIVA,LINAN,SCHME)Lophotrochozoa,((((((DROME,ANOGA),TRICA),DAPPU),STRMA),(IXOSC,PARTE)),CAEEL)Ecdysozoa)Protostomia,HOFMI)Bilateria,(((NEMVE,EXAPA),ACRMI),((HYDVU,CLYHE),AURAU))Cnidaria)Planulozoa,(HOIHO,TRIAD)Placozoa)Parahoxozoa,(PLEBA,MNELE)Ctenophora)Ctenophora_Parahoxozoa,(SYCCI,AMPQU)Porifera)Metazoa,SALRO)Choanozoa,CAPOW)Filozoa;
```
* XAMB
```
(((((((((((((((HOMSA,MUSMU),(CHEMY,GALGA)),XENTR),LATCH),(((MAYZE,HIPCO),DANRE),LEPOC)),CALMI)Vertebrata,CIOIN)Olfactores,BRALA)Chordata,(HOFMI,((SACKO,PTYFL),(STRPU,ACAPL))Ambulacraria)Xenambulacraria)Deuterostomia,(((EUPSC,LOTGI,(MIZYE,CRAGI))Mollusca,(CAPTE,HELRO),ADIVA,LINAN,SCHME)Lophotrochozoa,((((((DROME,ANOGA),TRICA),DAPPU),STRMA),(IXOSC,PARTE)),CAEEL)Ecdysozoa)Protostomia)Bilateria,(((NEMVE,EXAPA),ACRMI),((HYDVU,CLYHE),AURAU))Cnidaria)Planulozoa,(HOIHO,TRIAD)Placozoa)Parahoxozoa,(PLEBA,MNELE)Ctenophora,(SYCCI,AMPQU)Porifera)Metazoa,SALRO)Choanozoa,CAPOW)Filozoa;
```

* XNEP
```
((((((((((((((((HOMSA,MUSMU),(CHEMY,GALGA)),XENTR),LATCH),(((MAYZE,HIPCO),DANRE),LEPOC)),CALMI)Vertebrata,CIOIN)Olfactores,BRALA)Chordata,((SACKO,PTYFL),(STRPU,ACAPL))Ambulacraria)Deuterostomia,(((EUPSC,LOTGI,(MIZYE,CRAGI))Mollusca,(CAPTE,HELRO),ADIVA,LINAN,SCHME)Lophotrochozoa,((((((DROME,ANOGA),TRICA),DAPPU),STRMA),(IXOSC,PARTE)),CAEEL)Ecdysozoa)Protostomia)Nephrozoa,HOFMI)Bilateria,(((NEMVE,EXAPA),ACRMI),((HYDVU,CLYHE),AURAU))Cnidaria)Planulozoa,(HOIHO,TRIAD)Placozoa)Parahoxozoa,(PLEBA,MNELE)Ctenophora,(SYCCI,AMPQU)Porifera)Metazoa,SALRO)Choanozoa,CAPOW)Filozoa;
```

CONS is the tree we want to use for everything (since  alternate topologies generate similar results).
CS and PS will be used only for checking that the polytomies in our tree have no incidence on the estimation of metazoan novelties.
XNEP and XAMB are to be used for checking that the polytomies in our tree don't affect estimation of novelties in  Bilateria/Nephrozoa (LCA of protostomes and deuterostomes...) and deuterostomia LCA.
And ofc, the initial multifurcating tree, polytomic at the base of metazoa and at the base of bilateria

`filt.clusters` is filtered out version of the total multi_species block (minus cluster 1, 8406 blocks considered as one multi species block, due to fusing orthologous groups together. Randomized chromosomes don't have such artifacts, no need to filter them.

* Ctenophore-sister
```
cd 02_gene_density_analysis
python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/cteno_sister_bila_poly.tre -n Choanozoa Metazoa Porifera_Parahoxozoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel

#clusters list, to draw Venn diagrams of clusters found in nodes
python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/cteno_sister_bila_poly.tre -n Metazoa Porifera_Parahoxozoa -m 2 -r clusters_list -t novel|cut -f1,2 > venn_diag/CS.tsv



#counts for randomized blocks
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.1/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.1/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/cteno_sister_bila_poly.tre -n Choanozoa Metazoa Porifera_Parahoxozoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.2/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.2/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/cteno_sister_bila_poly.tre -n Choanozoa Metazoa Porifera_Parahoxozoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.3/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.3/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/cteno_sister_bila_poly.tre -n Choanozoa Metazoa Porifera_Parahoxozoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel

```

Ctenophore-sister    | ancestral | novel | recency
---------------------|-----------|-------|--------
Choanozoa            |    17     |   14  |    1
Metazoa              |     0     |    9  |    2
Porifera_Parahoxozoa |     9     |   25  |    3
Parahoxozoa          |    50     |    6  |    4
Planulozoa           |    63     |  162  |    5
Cnidaria             |   254     |   11  |    6
Bilateria            |   335     |  256  |    6
Deuterostomia        |   545     |    3  |    7
Protostomia          |   612     |   16  |    7
Chordata             |   511     |    1  |    8
Ambulacraria         |   430     |    6  |    8
Lophotrochozoa       |   798     |   91  |    8
Ecdysozoa            |   241     |    1  |    8
Olfactores           |   147     |    1  |    9
Mollusca             |   662     |   66  |    9
Vertebrata           |    70     |  170  |   10



* Porifera-sister
```
python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/sponge_sister_bila_poly.tre -n Choanozoa Metazoa Ctenophora_Parahoxozoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel

#clusters list, to draw Venn diagrams of clusters found in nodes
python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/sponge_sister_bila_poly.tre -n Metazoa Ctenophora_Parahoxozoa -m 2 -r clusters_list -t novel|cut -f1,2 > venn_diag/PS.tsv

#counts for randomized blocks
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.1/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.1/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/sponge_sister_bila_poly.tre -n Choanozoa Metazoa Ctenophora_Parahoxozoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.2/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.2/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/sponge_sister_bila_poly.tre -n Choanozoa Metazoa Ctenophora_Parahoxozoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.3/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.3/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/sponge_sister_bila_poly.tre -n Choanozoa Metazoa Ctenophora_Parahoxozoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
```

Porifera-sister        | ancestral | novel | recency
-----------------------|-----------|-------|--------
Choanozoa              |    17     |   14  |    1
Metazoa                |     0     |   25  |    2
Ctenophora_Parahoxozoa |    31     |    7  |    3
Parahoxozoa            |    50     |    6  |    4
Planulozoa             |    63     |  162  |    5
Cnidaria               |   254     |   11  |    6
Bilateria              |   335     |  256  |    6
Deuterostomia          |   545     |    3  |    7
Protostomia            |   612     |   16  |    7
Chordata               |   511     |    1  |    8
Ambulacraria           |   430     |    6  |    8
Lophotrochozoa         |   798     |   91  |    8
Ecdysozoa              |   241     |    1  |    8
Olfactores             |   147     |    1  |    9
Mollusca               |   662     |   66  |    9
Vertebrata             |    70     |  170  |   10



* Xenacoelomorpha

```
python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/meta_poly_xenacoelomorpha.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Nephrozoa Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel

#clusters list, to draw Venn diagrams of clusters found in nodes
python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/meta_poly_xenacoelomorpha.tre -n Bilateria Nephrozoa Deuterostomia Ambulacraria -m 2 -r clusters_list -t novel|cut -f1,2 > venn_diag/XNEP.tsv

#count blocks in randomized genomes
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.1/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.1/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/meta_poly_xenacoelomorpha.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Nephrozoa Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.2/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.2/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/meta_poly_xenacoelomorpha.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Nephrozoa Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.3/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.3/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/meta_poly_xenacoelomorpha.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Nephrozoa Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel

```

Xenacoelomorpha        | ancestral | novel | recency
-----------------------|-----------|-------|--------
Choanozoa              |    17     |   14  |    1
Metazoa                |     0     |   34  |    2
Parahoxozoa            |    50     |    6  |    3
Planulozoa             |    60     |  162  |    4
Cnidaria               |   254     |   11  |    5
Bilateria              |   335     |   38  |    5
Nephrozoa              |   354     |  224  |    6
Deuterostomia          |   545     |    3  |    7
Protostomia            |   612     |   16  |    7
Chordata               |   511     |    1  |    8
Ambulacraria           |   430     |    6  |    8
Lophotrochozoa         |   798     |   91  |    8
Ecdysozoa              |   241     |    1  |    8
Olfactores             |   147     |    1  |    9
Mollusca               |   662     |   66  |    9
Vertebrata             |    70     |  170  |   10

* Xenambulacraria

```
python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/meta_poly_xenambulacraria.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Xenambulacraria Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel

#clusters list, to draw Venn diagrams of clusters found in nodes
python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/meta_poly_xenambulacraria.tre -n Bilateria Deuterostomia Xenambulacraria Ambulacraria -m 2 -r clusters_list -t novel|cut -f1,2 > venn_diag/XAMB.tsv

#Isolate blocks, and clusters



#random block counts
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.1/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.1/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/meta_poly_xenambulacraria.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Xenambulacraria Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.2/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.2/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/meta_poly_xenambulacraria.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Xenambulacraria Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.3/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.3/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/meta_poly_xenambulacraria.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Xenambulacraria Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel


```
No syntenic novelty retained in xenacoelomorpha

Xenambulacraria        | ancestral | novel | recency
-----------------------|-----------|-------|--------
Choanozoa              |    17     |   14  |    1
Metazoa                |     0     |   34  |    2
Parahoxozoa            |    50     |    6  |    3
Planulozoa             |   63     |  162  |    4
Cnidaria               |   254     |   11  |    5
Bilateria              |   328     |  237  |    5
Deuterostomia          |   568     |    3  |    6
Protostomia            |   612     |   16  |    6
Xenambulacraria        |   511     |    0  |    7
Chordata               |   517     |    1  |    8
Ambulacraria           |   430     |    6  |    8
Lophotrochozoa         |   798     |   91  |    8
Ecdysozoa              |   241     |    1  |    8
Olfactores             |   147     |    1  |    9
Mollusca               |   662     |   66  |    9
Vertebrata             |    70     |  170  |   10


* Consensus
```
python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/meta_poly_bila_poly.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel

#clusters list, to draw Venn diagrams of clusters found in nodes
python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/meta_poly_bila_poly.tre -n Metazoa Bilateria Deuterostomia Ambulacraria -m 2 -r clusters_list -t novel|cut -f1,2 > venn_diag/CONS.tsv


nodelist=(Bilateria Lophotrochozoa Metazoa Parahoxozoa Planulozoa Vertebrata)
#Block types
for node in "${nodelist[@]}"; do
    python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/meta_poly_bila_poly.tre -n $node -m 2 -r blocks_list -t novel|cut -f2- > bynode/CONS/$node.novel.synt
    python3 BlocksByNode.py -c filt.clusters -b ../01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -s trees/meta_poly_bila_poly.tre -n $node -m 2 -r clusters_list -t novel|cut -f1,3- > bynode/CONS/$node.novel.clusters;
    done

#random block counts
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.1/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.1/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/meta_poly_bila_poly.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.2/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.2/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/meta_poly_bila_poly.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
python3 BlocksByNode.py -c ../01_microsynteny/randomized_chrom/rand.3/Orthofinder/5.blocks.3.syn.clusters -b ../01_microsynteny/randomized_chrom/rand.3/Orthofinder/5.blocks.3.syn_corrected.synt -s trees/meta_poly_bila_poly.tre -n Choanozoa Metazoa Parahoxozoa Planulozoa Cnidaria Bilateria Deuterostomia Protostomia Chordata Ambulacraria Lophotrochozoa Ecdysozoa Olfactores Mollusca Vertebrata -m 2 -r short -t ancestral novel
```

Consensus              | ancestral | novel | recency
-----------------------|-----------|-------|--------
Choanozoa              |    17     |   14  |    1
Metazoa                |     0     |   34  |    2
Parahoxozoa            |    50     |    6  |    3
Planulozoa             |    63     |  162  |    4
Cnidaria               |   254     |   11  |    5
Bilateria              |   335     |  256  |    5
Deuterostomia          |   545     |    3  |    6
Protostomia            |   612     |   16  |    6
Chordata               |   511     |    1  |    7
Ambulacraria           |   430     |    6  |    7
Lophotrochozoa         |   798     |   91  |    7
Ecdysozoa              |   241     |    1  |    7
Olfactores             |   147     |    1  |    8
Mollusca               |   662     |   66  |    8
Vertebrata             |    70     |  170  |    9




### 5.2 Measuring gene density at key nodes
We'll use consensus tree

```
cd 02_gene_density_analysis/bynode/CONS
#randomized blocks
nodelist=(Bilateria Lophotrochozoa Metazoa Parahoxozoa Planulozoa Vertebrata)
for node in "${nodelist[@]}"; do
    while read species; do
    echo "Rscript pick_random_blocks.R --only-genome=$species -n 100 random/$node/ $node.novel.synt ../../../01_microsynteny/chrom/$species.chrom";
    done < $node.species;
done > random.sh

cat random.sh |slurmtasks -n random_blocks |sbatch

for node in "${nodelist[@]}"; do
    cat random/$node/*tsv|grep -v "block_id" > $node.novel.random.synt ;
done
```

Now for measuring gene density in every animal and node, we use `make_tidy_density_df`.  It can take a while (2-3 hours for the 49 species dataset, 6 different syntenic nodes, everything runs sequentially on a single core of an AMD Opteron 6320).
This outputs a huge table with gene density by block, genome size, median intergenic distance, node, taxon, species, gene list, para/not_para classification (>40% or less than  40% of genes belonging to the same OG) (i.e. all info needed to build figures).

Total assembly size is measured by parsing fasta files
```
python3 make_tidy_density_df.py -s ../bynode/CONS/ -g ../genomes/ -c ../../01_microsynteny/chrom/ -m ../../01_microsynteny/chrom_of/5.blocks.3.syn.clusters -og ../../01_microsynteny/Orthofinder.clus
```
Make data to be used for scatterplots of SF4 (run in the same folder as tidy density df)
```
figure2_prep_data_scatterplots.py
```## 6 GO_terms enrichment analysis
We'll use eggnogmapper.
Two steps: first emapper is homology searches with diamond, second is the actual annotation
```
fastafolder="00_orthology_assignment/proteins_processed/"
outfolder="05_GO/emapper"

cd $fastafolder
for f in *.fasta; do echo "python /apps/eggnogmapper/2.0.0/emapper.py -m diamond --no_annot --no_file_comments --cpu 1 -i 00_orthology_assignment/proteins_processed/$f -o 05_GO/emapper/$f"; done>$outfolder/homology_searches
cd $outfolder

cat homology_searches |slurmtasks -n emapper -m 10 |sbatch

ls |grep emapper.seed_orthologs|cut -f1 -d '.' > species.list
while read species; do echo "emapper.py --annotate_hits_table $species.fasta.emapper.seed_orthologs --no_file_comments -o $species --cpu 2"; done<species_list>annotation_job
cat annotation_job |slurmtasks -n emapper -m 1 -c 2 -f array-2core |sbatch

```
We want to make ids2go file for  all our species. We'll change the comma-delimited GO IDs to semi-colon delimited, as this is the format the GOATOOLS API parses.
```
for species in $(cat species.list); do cut -f 1,7 $species.emapper.annotations|grep GO:|sed 's/,/;/g' > ids2GO_files/$species.ids2go; done
```
Isolate lists of proteins for each species
First, we make lists for all the species of all the proteins in each proteome
```
cd 00_orthology_assignment/proteins_processed/
species_list=$(ls *.fasta|cut -f1 -d '.')
for species in $species_list; do
    grep '>' $species.fasta | cut -f 2 -d '>' > 05_GO/emapper/ids2GO_files/$species.list
done

```
we won't do parahoxozoa enrichments, since blocks we found are roughly as many as background ones 
```

cd 05_GO/GO_enrichment

annotations="../emapper/annotations"
species_list=$(ls $annotations/*emapper.annotations|rev|cut -d '/' -f 1|rev|cut -f1 -d '.') #cut to keep only the last field
ids2go_folder=../emapper/ids2GO_files

#MLCA enrichments
name=Metazoa
syntfile=../../02_REDUX_gene_density_analysis/bynode/CONS/Metazoa.novel.synt
for species in $species_list; do
    grep -P "^\d+\t$species" $syntfile | cut -f10|tr "," "\n"> Metazoa/$species.$name.list;
done

#PLCA enrichments
name=Planulozoa
syntfile=../../02_REDUX_gene_density_analysis/bynode/CONS/Planulozoa.novel.synt
for species in $species_list; do
    grep -P "^\d+\t$species" $syntfile | cut -f10|tr "," "\n"> Planulozoa/$species.$name.list;
done

#BLCA enrichments
name=Bilateria
syntfile=../../02_REDUX_gene_density_analysis/bynode/CONS/Bilateria.novel.synt
for species in $species_list; do
    grep -P "^\d+\t$species" $syntfile | cut -f10|tr "," "\n"> Bilateria/$species.$name.list;
done

#remove empty files (no blocks)
find */* -empty -type f -delete

#get latest version of GO terms (downloaded on the 5th april 2020)
wget http://purl.obolibrary.org/obo/go/go-basic.obo

#run GO enrichment for all species posessing MCLA blocks
species_list=$(ls Metazoa |cut -f2 -d "/"|cut -f1 -d ".") #only species in given folder
for species in $species_list; do
    ids2go=$ids2go_folder/$species.ids2go
    background_seq=$ids2go_folder/$species.list
    sample_seq=Metazoa/$species.Metazoa.list
    python3 /proj/robert/scripts/GO/GO_analysis.py -i $ids2go -go go-basic.obo -b $background_seq -s $sample_seq;
done

#run GO enrichment for all species posessing PCLA blocks
species_list=$(ls Planulozoa |cut -f2 -d "/"|cut -f1 -d ".") #only species in given folder
for species in $species_list; do
    ids2go=$ids2go_folder/$species.ids2go
    background_seq=$ids2go_folder/$species.list
    sample_seq=Planulozoa/$species.Planulozoa.list
    python3 /proj/robert/scripts/GO/GO_analysis.py -i $ids2go -go go-basic.obo -b $background_seq -s $sample_seq;
done

#run GO enrichment for all species posessing BLCA blocks
species_list=$(ls Bilateria |cut -f2 -d "/"|cut -f1 -d ".") #only species in given folder
for species in $species_list; do
    ids2go=$ids2go_folder/$species.ids2go
    background_seq=$ids2go_folder/$species.list
    sample_seq=Bilateria/$species.Bilateria.list
    python3 /proj/robert/scripts/GO/GO_analysis.py -i $ids2go -go go-basic.obo -b $background_seq -s $sample_seq;
done
```

check which terms are enriched in which taxons. GO enrichment results are in the GO_enrichment_results folder
```
cd 05_GO/GO_enrichment
enrichment_comparisons.py
```## 7. Block correlation
### 7.1 Download Expression data for CALMI and MUSMU
All other transcript abundances taken from Zieger et al. 2020 (doi: 10.1016/j.cub.2020.10.004)
#### CALMI
We first want to use our filtered gff file (no MT genome). And download the transcripts
```
08_exp/CALMI/CALMI.gff

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/165/045/GCF_000165045.1_Callorhinchus_milii-6.1.3/GCF_000165045.1_Callorhinchus_milii-6.1.3_rna.fna.gz -O CALMI.transcripts.fa.gz
```

Now this is for renaming the transcripts in the mRNA file by the protein names and filter out transcripts that are not the ones encoding our filtered proteins
```
/proj/robert/scripts/dev_RNAseq/rename_filter_transcripts.py -a CALMI_NCBI.gff -ft mRNA -fp CDS -kpp Parent -kpt ID -kp protein_id -kt transcript_id -o CALMI_protnames.fna -f CALMI.transcripts.fa -c 01_microsynteny/chrom/CALMI.chrom


```

Download all the reads
```
module load sratoolkit
fasterq-dump SRR1735385 -e 8 -o thymus
fasterq-dump SRR513757 -e 8 -o testis
fasterq-dump SRR513758 -e 8 -o spleen
fasterq-dump SRR513759 -e 8 -o ovary
fasterq-dump SRR513760 -e 8 -o liver
fasterq-dump SRR514104 -e 8 -o muscle
fasterq-dump SRR514105 -e 8 -o kidney
fasterq-dump SRR514106 -e 8 -o intestine
fasterq-dump SRR514107 -e 8 -o heart
fasterq-dump SRR514109 -e 8 -o brain
fasterq-dump SRR534176 -e 8 -o gills
```
Now for the kallisto run
```
module load kallisto
#cerate index
kallisto index -i CALMI.index CALMI_protnames.fna

reads=08_exp/CALMI/reads

le1="_1.fastq"
le2="_2.fastq"

#submit kallisto quantification jobs
for file in $(ls reads/|rev|cut -f 2- -d '_'|rev|uniq); do
    echo "kallisto quant --index=CALMI.index --rf-stranded --output-dir=$file --plaintext reads/$file$le1  reads/$file$le2" ;
done |
/proj/rpz/slurm_scripts/slurmtasks --mem 5 --name quantCALMI |
sbatch
```


```
paste */abundance.tsv |  cut -f 1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115 > CALMI_transcript_tpms_all_samples.tsv
ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id$_\n"' > header.tsv
cat header.tsv CALMI_transcript_tpms_all_samples.tsv | grep -v "tpm" > CALMI_transcript_tpms_all_samples.tsv2
mv CALMI_transcript_tpms_all_samples.tsv2 CALMI_transcript_tpms_all_samples.tsv
rm -f header.tsv
```


## MUSMU
```
python3 /proj/robert/scripts/dev_RNAseq/rename_filter_transcripts.py -a 00_orthology_assignment/coordinates/gff/MUSMU_filt_NCBI.gff3 -ft mRNA -fp CDS -kpp gene -kpt gene -kp protein_id -kt Name -o MUSMU_protnames.fna -f MUSMU.transcripts.fa -c 01_microsynteny/chrom/MUSMU.chrom
```

```
fasterq-dump SRR5273648 -e 4 -o female_bone_marrow_a
fasterq-dump SRR5273664 -e 4 -o female_bone_marrow_b


module load sratoolkit
prefetch SRR5273654
cd SRR5273654
fastq-dump --split-e --skip-technical SRR5273654.sra
mv SRR5273654_1.fastq ../female_adrenal_gland_b_1.fastq
mv SRR5273654_2.fastq ../female_adrenal_gland_b_2.fastq

prefetch SRR5273670
cd SRR5273670
fastq-dump --split-e --skip-technical SRR5273670.sra
mv SRR5273670_1.fastq ../female_adrenal_gland_a_1.fastq
mv SRR5273670_2.fastq ../female_adrenal_gland_a_2.fastq
cd ../
rm -r SRR5273670

prefetch SRR5273635
cd SRR5273635
fastq-dump --split-e --skip-technical SRR5273635.sra
mv SRR5273635_1.fastq ../female_brain_a_1.fastq
mv SRR5273635_2.fastq ../female_brain_a_2.fastq
cd ../
rm -r SRR5273635

prefetch SRR5273637
cd SRR5273637
fastq-dump --split-e --skip-technical SRR5273637.sra
mv SRR5273637_1.fastq ../female_brain_b_1.fastq
mv SRR5273637_2.fastq ../female_brain_b_2.fastq
cd ../
rm -r SRR5273637

#todo
prefetch SRR5273657
cd SRR5273657
fastq-dump --split-e --skip-technical SRR5273657.sra
mv SRR5273657_1.fastq ../female_brain_c_1.fastq
mv SRR5273657_2.fastq ../female_brain_c_2.fastq
cd ../
rm -r SRR5273657

prefetch SRR5273673
cd SRR5273673
fastq-dump --split-e --skip-technical SRR5273673.sra
mv SRR5273673_1.fastq ../female_brain_d_1.fastq
mv SRR5273673_2.fastq ../female_brain_d_2.fastq
cd ../
rm -r SRR5273673


prefetch SRR5273646
cd SRR5273646
fastq-dump --split-e --skip-technical SRR5273646.sra
mv SRR5273646_1.fastq ../female_forestomach_a_1.fastq
mv SRR5273646_2.fastq ../female_forestomach_a_2.fastq
cd ../
rm -r SRR5273646


prefetch SRR5273662
cd SRR5273662
fastq-dump --split-e --skip-technical SRR5273662.sra
mv SRR5273662_1.fastq ../female_forestomach_b_1.fastq
mv SRR5273662_2.fastq ../female_forestomach_b_2.fastq
cd ../
rm -r SRR5273662

prefetch SRR5273651
cd SRR5273651
fastq-dump --split-e --skip-technical SRR5273651.sra
mv SRR5273651_1.fastq ../female_heart_a_1.fastq
mv SRR5273651_2.fastq ../female_heart_a_2.fastq
cd ../
rm -r SRR5273651

prefetch SRR5273667
cd SRR5273667
fastq-dump --split-e --skip-technical SRR5273667.sra
mv SRR5273667_1.fastq ../female_heart_b_1.fastq
mv SRR5273667_2.fastq ../female_heart_b_2.fastq
cd ../
rm -r SRR5273667


prefetch SRR5273655
cd SRR5273655
fastq-dump --split-e --skip-technical SRR5273655.sra
mv SRR5273655_1.fastq ../female_kidney_a_1.fastq
mv SRR5273655_2.fastq ../female_kidney_a_2.fastq
cd ../
rm -r SRR5273655

prefetch SRR5273671
cd SRR5273671
fastq-dump --split-e --skip-technical SRR5273671.sra
mv SRR5273671_1.fastq ../female_kidney_b_1.fastq
mv SRR5273671_2.fastq ../female_kidney_b_2.fastq
cd ../
rm -r SRR5273671

prefetch SRR5273644
cd SRR5273644
fastq-dump --split-e --skip-technical SRR5273644.sra
mv SRR5273644_1.fastq ../female_large_intestine_a_1.fastq
mv SRR5273644_2.fastq ../female_large_intestine_a_2.fastq
cd ../
rm -r SRR5273644

prefetch SRR5273660
cd SRR5273660
fastq-dump --split-e --skip-technical SRR5273660.sra
mv SRR5273660_1.fastq ../female_large_intestine_b_1.fastq
mv SRR5273660_2.fastq ../female_large_intestine_b_2.fastq
cd ../
rm -r SRR5273660

prefetch SRR5273634
cd SRR5273634
fastq-dump --split-e --skip-technical SRR5273634.sra
mv SRR5273634_1.fastq ../female_liver_a_1.fastq
mv SRR5273634_2.fastq ../female_liver_a_2.fastq
cd ../
rm -r SRR5273634

prefetch SRR5273636
cd SRR5273636
fastq-dump --split-e --skip-technical SRR5273636.sra
mv SRR5273636_1.fastq ../female_liver_b_1.fastq
mv SRR5273636_2.fastq ../female_liver_b_2.fastq
cd ../
rm -r SRR5273636

prefetch SRR5273656
cd SRR5273656
fastq-dump --split-e --skip-technical SRR5273656.sra
mv SRR5273656_1.fastq ../female_liver_c_1.fastq
mv SRR5273656_2.fastq ../female_liver_c_2.fastq
cd ../
rm -r SRR5273656

prefetch SRR5273672
cd SRR5273672
fastq-dump --split-e --skip-technical SRR5273672.sra
mv SRR5273672_1.fastq ../female_liver_d_1.fastq
mv SRR5273672_2.fastq ../female_liver_d_2.fastq
cd ../
rm -r SRR5273672

prefetch SRR5273652
cd SRR5273652
fastq-dump --split-e --skip-technical SRR5273652.sra
mv SRR5273652_1.fastq ../female_lung_a_1.fastq
mv SRR5273652_2.fastq ../female_lung_a_2.fastq
cd ../
rm -r SRR5273652

prefetch SRR5273668
cd SRR5273668
fastq-dump --split-e --skip-technical SRR5273668.sra
mv SRR5273668_1.fastq ../female_lung_b_1.fastq
mv SRR5273668_2.fastq ../female_lung_b_2.fastq
cd ../
rm -r SRR5273668

prefetch SRR5273643
cd SRR5273643
fastq-dump --split-e --skip-technical SRR5273643.sra
mv SRR5273643_1.fastq ../female_muscle_a_1.fastq
mv SRR5273643_2.fastq ../female_muscle_a_2.fastq
cd ../
rm -r SRR5273643

prefetch SRR5273659
cd SRR5273659
fastq-dump --split-e --skip-technical SRR5273659.sra
mv SRR5273659_1.fastq ../female_muscle_b_1.fastq
mv SRR5273659_2.fastq ../female_muscle_b_2.fastq
cd ../
rm -r SRR5273659

prefetch SRR5273649
cd SRR5273649
fastq-dump --split-e --skip-technical SRR5273649.sra
mv SRR5273649_1.fastq ../female_ovary_a_1.fastq
mv SRR5273649_2.fastq ../female_ovary_a_2.fastq
cd ../
rm -r SRR5273649

prefetch SRR5273665
cd SRR5273665
fastq-dump --split-e --skip-technical SRR5273665.sra
mv SRR5273665_1.fastq ../female_ovary_b_1.fastq
mv SRR5273665_2.fastq ../female_ovary_b_2.fastq
cd ../
rm -r SRR5273665

prefetch SRR5273645
cd SRR5273645
fastq-dump --split-e --skip-technical SRR5273645.sra
mv SRR5273645_1.fastq ../female_small_intestine_a_1.fastq
mv SRR5273645_2.fastq ../female_small_intestine_a_2.fastq
cd ../
rm -r SRR5273645

prefetch SRR5273661
cd SRR5273661
fastq-dump --split-e --skip-technical SRR5273661.sra
mv SRR5273661_1.fastq ../female_small_intestine_b_1.fastq
mv SRR5273661_2.fastq ../female_small_intestine_b_2.fastq
cd ../
rm -r SRR5273661

prefetch SRR5273653
cd SRR5273653
fastq-dump --split-e --skip-technical SRR5273653.sra
mv SRR5273653_1.fastq ../female_spleen_a_1.fastq
mv SRR5273653_2.fastq ../female_spleen_a_2.fastq
cd ../
rm -r SRR5273653

prefetch SRR5273669
cd SRR5273669
fastq-dump --split-e --skip-technical SRR5273669.sra
mv SRR5273669_1.fastq ../female_spleen_b_1.fastq
mv SRR5273669_2.fastq ../female_spleen_b_2.fastq
cd ../
rm -r SRR5273669

prefetch SRR5273647
cd SRR5273647
fastq-dump --split-e --skip-technical SRR5273647.sra
mv SRR5273647_1.fastq ../female_stomach_a_1.fastq
mv SRR5273647_2.fastq ../female_stomach_a_2.fastq
cd ../
rm -r SRR5273647

prefetch SRR5273663
cd SRR5273663
fastq-dump --split-e --skip-technical SRR5273663.sra
mv SRR5273663_1.fastq ../female_stomach_b_1.fastq
mv SRR5273663_2.fastq ../female_stomach_b_2.fastq
cd ../
rm -r SRR5273663

prefetch SRR5273650
cd SRR5273650
fastq-dump --split-e --skip-technical SRR5273650.sra
mv SRR5273650_1.fastq ../female_thymus_a_1.fastq
mv SRR5273650_2.fastq ../female_thymus_a_2.fastq
cd ../
rm -r SRR5273650

prefetch SRR5273666
cd SRR5273666
fastq-dump --split-e --skip-technical SRR5273666.sra
mv SRR5273666_1.fastq ../female_thymus_b_1.fastq
mv SRR5273666_2.fastq ../female_thymus_b_2.fastq
cd ../
rm -r SRR5273666

prefetch SRR5273642
cd SRR5273642
fastq-dump --split-e --skip-technical SRR5273642.sra
mv SRR5273642_1.fastq ../female_uterus_1.fastq
mv SRR5273642_2.fastq ../female_uterus_2.fastq
cd ../
rm -r SRR5273642

prefetch SRR5273658
cd SRR5273658
fastq-dump --split-e --skip-technical SRR5273658.sra
mv SRR5273658_1.fastq ../female_vesicular_gland_1.fastq
mv SRR5273658_2.fastq ../female_vesicular_gland_2.fastq
cd ../
rm -r SRR5273658

prefetch SRR5273686
cd SRR5273686
fastq-dump --split-e --skip-technical SRR5273686.sra
mv SRR5273686_1.fastq ../male_adrenal_gland_a_1.fastq
mv SRR5273686_2.fastq ../male_adrenal_gland_a_2.fastq
cd ../
rm -r SRR5273686

prefetch SRR5273702
cd SRR5273702
fastq-dump --split-e --skip-technical SRR5273702.sra
mv SRR5273702_1.fastq ../male_adrenal_gland_b_1.fastq
mv SRR5273702_2.fastq ../male_adrenal_gland_b_2.fastq
cd ../
rm -r SRR5273702

prefetch SRR5273680
cd SRR5273680
fastq-dump --split-e --skip-technical SRR5273680.sra
mv SRR5273680_1.fastq ../male_bone_marrow_a_1.fastq
mv SRR5273680_2.fastq ../male_bone_marrow_a_2.fastq
cd ../
rm -r SRR5273680

prefetch SRR5273696
cd SRR5273696
fastq-dump --split-e --skip-technical SRR5273696.sra
mv SRR5273696_1.fastq ../male_bone_marrow_b_1.fastq
mv SRR5273696_2.fastq ../male_bone_marrow_b_2.fastq
cd ../
rm -r SRR5273696

prefetch SRR5273639
cd SRR5273639
fastq-dump --split-e --skip-technical SRR5273639.sra
mv SRR5273639_1.fastq ../male_brain_a_1.fastq
mv SRR5273639_2.fastq ../male_brain_a_2.fastq
cd ../
rm -r SRR5273639

prefetch SRR5273641
cd SRR5273641
fastq-dump --split-e --skip-technical SRR5273641.sra
mv SRR5273641_1.fastq ../male_brain_b_1.fastq
mv SRR5273641_2.fastq ../male_brain_b_2.fastq
cd ../
rm -r SRR5273641

prefetch SRR5273689
cd SRR5273689
fastq-dump --split-e --skip-technical SRR5273689.sra
mv SRR5273689_1.fastq ../male_brain_c_1.fastq
mv SRR5273689_2.fastq ../male_brain_c_2.fastq
cd ../
rm -r SRR5273689

prefetch SRR5273705
cd SRR5273705
fastq-dump --split-e --skip-technical SRR5273705.sra
mv SRR5273705_1.fastq ../male_brain_d_1.fastq
mv SRR5273705_2.fastq ../male_brain_d_2.fastq
cd ../
rm -r SRR5273705

prefetch SRR5273678
cd SRR5273678
fastq-dump --split-e --skip-technical SRR5273678.sra
mv SRR5273678_1.fastq ../male_forestomach_a_1.fastq
mv SRR5273678_2.fastq ../male_forestomach_a_2.fastq
cd ../
rm -r SRR5273678

prefetch SRR5273694
cd SRR5273694
fastq-dump --split-e --skip-technical SRR5273694.sra
mv SRR5273694_1.fastq ../male_forestomach_b_1.fastq
mv SRR5273694_2.fastq ../male_forestomach_b_2.fastq
cd ../
rm -r SRR5273694

prefetch SRR5273683
cd SRR5273683
fastq-dump --split-e --skip-technical SRR5273683.sra
mv SRR5273683_1.fastq ../male_heart_a_1.fastq
mv SRR5273683_2.fastq ../male_heart_a_2.fastq
cd ../
rm -r SRR5273683

prefetch SRR5273699
cd SRR5273699
fastq-dump --split-e --skip-technical SRR5273699.sra
mv SRR5273699_1.fastq ../male_heart_b_1.fastq
mv SRR5273699_2.fastq ../male_heart_b_2.fastq
cd ../
rm -r SRR5273699

prefetch SRR5273687
cd SRR5273687
fastq-dump --split-e --skip-technical SRR5273687.sra
mv SRR5273687_1.fastq ../male_kidney_a_1.fastq
mv SRR5273687_2.fastq ../male_kidney_a_2.fastq
cd ../
rm -r SRR5273687

prefetch SRR5273703
cd SRR5273703
fastq-dump --split-e --skip-technical SRR5273703.sra
mv SRR5273703_1.fastq ../male_kidney_b_1.fastq
mv SRR5273703_2.fastq ../male_kidney_b_2.fastq
cd ../
rm -r SRR5273703

prefetch SRR5273676
cd SRR5273676
fastq-dump --split-e --skip-technical SRR5273676.sra
mv SRR5273676_1.fastq ../male_large_intestine_a_1.fastq
mv SRR5273676_2.fastq ../male_large_intestine_a_2.fastq
cd ../
rm -r SRR5273676

prefetch SRR5273692
cd SRR5273692
fastq-dump --split-e --skip-technical SRR5273692.sra
mv SRR5273692_1.fastq ../male_large_intestine_b_1.fastq
mv SRR5273692_2.fastq ../male_large_intestine_b_2.fastq
cd ../
rm -r SRR5273692

prefetch SRR5273638
cd SRR5273638
fastq-dump --split-e --skip-technical SRR5273638.sra
mv SRR5273638_1.fastq ../male_liver_a_1.fastq
mv SRR5273638_2.fastq ../male_liver_a_2.fastq
cd ../
rm -r SRR5273638

prefetch SRR5273640
cd SRR5273640
fastq-dump --split-e --skip-technical SRR5273640.sra
mv SRR5273640_1.fastq ../male_liver_b_1.fastq
mv SRR5273640_2.fastq ../male_liver_b_2.fastq
cd ../
rm -r SRR5273640

prefetch SRR5273688
cd SRR5273688
fastq-dump --split-e --skip-technical SRR5273688.sra
mv SRR5273688_1.fastq ../male_liver_c_1.fastq
mv SRR5273688_2.fastq ../male_liver_c_2.fastq
cd ../
rm -r SRR5273688

prefetch SRR5273704
cd SRR5273704
fastq-dump --split-e --skip-technical SRR5273704.sra
mv SRR5273704_1.fastq ../male_liver_d_1.fastq
mv SRR5273704_2.fastq ../male_liver_d_2.fastq
cd ../
rm -r SRR5273704

prefetch SRR5273684
cd SRR5273684
fastq-dump --split-e --skip-technical SRR5273684.sra
mv SRR5273684_1.fastq ../male_lung_a_1.fastq
mv SRR5273684_2.fastq ../male_lung_a_2.fastq
cd ../
rm -r SRR5273684

prefetch SRR5273700
cd SRR5273700
fastq-dump --split-e --skip-technical SRR5273700.sra
mv SRR5273700_1.fastq ../male_lung_b_1.fastq
mv SRR5273700_2.fastq ../male_lung_b_2.fastq
cd ../
rm -r SRR5273700

prefetch SRR5273675
cd SRR5273675
fastq-dump --split-e --skip-technical SRR5273675.sra
mv SRR5273675_1.fastq ../male_muscle_a_1.fastq
mv SRR5273675_2.fastq ../male_muscle_a_2.fastq
cd ../
rm -r SRR5273675

prefetch SRR5273691
cd SRR5273691
fastq-dump --split-e --skip-technical SRR5273691.sra
mv SRR5273691_1.fastq ../male_muscle_b_1.fastq
mv SRR5273691_2.fastq ../male_muscle_b_2.fastq
cd ../
rm -r SRR5273691

prefetch SRR5273677
cd SRR5273677
fastq-dump --split-e --skip-technical SRR5273677.sra
mv SRR5273677_1.fastq ../male_small_intestine_a_1.fastq
mv SRR5273677_2.fastq ../male_small_intestine_a_2.fastq
cd ../
rm -r SRR5273677

prefetch SRR5273693
cd SRR5273693
fastq-dump --split-e --skip-technical SRR5273693.sra
mv SRR5273693_1.fastq ../male_small_intestine_b_1.fastq
mv SRR5273693_2.fastq ../male_small_intestine_b_2.fastq
cd ../
rm -r SRR5273693

prefetch SRR5273685
cd SRR5273685
fastq-dump --split-e --skip-technical SRR5273685.sra
mv SRR5273685_1.fastq ../male_spleen_a_1.fastq
mv SRR5273685_2.fastq ../male_spleen_a_2.fastq
cd ../
rm -r SRR5273685

prefetch SRR5273701
cd SRR5273701
fastq-dump --split-e --skip-technical SRR5273701.sra
mv SRR5273701_1.fastq ../male_spleen_b_1.fastq
mv SRR5273701_2.fastq ../male_spleen_b_2.fastq
cd ../
rm -r SRR5273701

prefetch SRR5273679
cd SRR5273679
fastq-dump --split-e --skip-technical SRR5273679.sra
mv SRR5273679_1.fastq ../male_stomach_a_1.fastq
mv SRR5273679_2.fastq ../male_stomach_a_2.fastq
cd ../
rm -r SRR5273679

prefetch SRR5273695
cd SRR5273695
fastq-dump --split-e --skip-technical SRR5273695.sra
mv SRR5273695_1.fastq ../male_stomach_b_1.fastq
mv SRR5273695_2.fastq ../male_stomach_b_2.fastq
cd ../
rm -r SRR5273695

prefetch SRR5273681
cd SRR5273681
fastq-dump --split-e --skip-technical SRR5273681.sra
mv SRR5273681_1.fastq ../male_testis_a_1.fastq
mv SRR5273681_2.fastq ../male_testis_a_2.fastq
cd ../
rm -r SRR5273681

prefetch SRR5273697
cd SRR5273697
fastq-dump --split-e --skip-technical SRR5273697.sra
mv SRR5273697_1.fastq ../male_testis_b_1.fastq
mv SRR5273697_2.fastq ../male_testis_b_2.fastq
cd ../
rm -r SRR5273697

prefetch SRR5273682
cd SRR5273682
fastq-dump --split-e --skip-technical SRR5273682.sra
mv SRR5273682_1.fastq ../male_thymus_a_1.fastq
mv SRR5273682_2.fastq ../male_thymus_a_2.fastq
cd ../
rm -r SRR5273682

prefetch SRR5273698
cd SRR5273698
fastq-dump --split-e --skip-technical SRR5273698.sra
mv SRR5273698_1.fastq ../male_thymus_b_1.fastq
mv SRR5273698_2.fastq ../male_thymus_b_2.fastq
cd ../
rm -r SRR5273698

prefetch SRR5273674
cd SRR5273674
fastq-dump --split-e --skip-technical SRR5273674.sra
mv SRR5273674_1.fastq ../male_vesicular_gland_a_1.fastq
mv SRR5273674_2.fastq ../male_vesicular_gland_a_2.fastq
cd ../
rm -r SRR5273674

prefetch SRR5273690
cd SRR5273690
fastq-dump --split-e --skip-technical SRR5273690.sra
mv SRR5273690_1.fastq ../male_vesicular_gland_b_1.fastq
mv SRR5273690_2.fastq ../male_vesicular_gland_b_2.fastq
cd ../
rm -r SRR5273690
```
kallisto run
```
module load kallisto
#cerate index
kallisto index -i MUSMU.index MUSMU_protnames.fna

reads=08_exp/MUSMU/reads

le1="_1.fastq"
le2="_2.fastq"

#submit kallisto quantification jobs
for file in $(ls reads/|rev|cut -f 2- -d '_'|rev|uniq); do
    echo "kallisto quant --index=MUSMU.index --rf-stranded --output-dir=$file --plaintext reads/$file$le1  reads/$file$le2" ;
done |
/proj/rpz/slurm_scripts/slurmtasks --mem 5 --name quantMUSMU |
sbatch



paste */abundance.tsv |  cut -f 1,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200,205,210,215,220,225,230,235,240,245,250,255,260,265,270,275,280,285,290,295,300,305,310,315,320,325,330,335,340,345,350,355,360,365,370,375,380,385,390 > MUSMU_transcript_tpms_all_samples.tsv
ls -1 */abundance.tsv | perl -ne 'chomp $_; if ($_ =~ /(\S+)\/abundance\.tsv/){print "\t$1"}' | perl -ne 'print "target_id$_\n"' > header.tsv
cat header.tsv MUSMU_transcript_tpms_all_samples.tsv | grep -v "tpm" > MUSMU_transcript_tpms_all_samples.tsv2
mv MUSMU_transcript_tpms_all_samples.tsv2 MUSMU_transcript_tpms_all_samples.tsv
rm -f header.tsv

python3 /proj/robert/scripts/Utilities/make_tpm_stage_medians.py "female_adrenal_gland female_bone_marrow female_brain female_forestomach female_heart female_kidney female_large_intestine female_liver female_lung female_muscle female_ovary female_small_intestine female_spleen female_stomach female_thymus female_uterus female_vesicular_gland male_adrenal_gland male_bone_marrow male_brain male_forestomach male_heart male_kidney male_large_intestine male_liver male_lung male_muscle male_small_intestine male_spleen male_stomach male_testis male_thymus male_vesicular_gland" -i MUSMU_transcript_tpms_all_samples.tsv

```

## Calculate the actual block correlations
Really large blocks in vertebrates. since it involve a bit of combinations, it can take a while, especially with the long gnathostome blocks.

Elephant shark
```
python3 /proj/robert/scripts/density/block_correlation_analysis.py -t 02_gene_density_analysis/density_whole_genome/key_nodes.tidydf.csv -s 01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -e CALMI/CALMI_transcript_tpms_all_samples.tsv -p CALMI -o CALMI_tidy_corr.csv 
```
Scallop
```
python3 /proj/robert/scripts/density/block_correlation_analysis.py -t 02_gene_density_analysis/density_whole_genome/key_nodes.tidydf.csv -s 01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -e MIZYE/MIZYE_transcript_tpms_all_samples_medians.tsv -p MIZYE -o MIZYE_tidy_corr.csv 
```
Mouse
```
python3 /proj/robert/scripts/density/block_correlation_analysis.py -t 02_gene_density_analysis/density_whole_genome/key_nodes.tidydf.csv -s 01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -e MUSMU/MUSMU_transcript_tpms_all_samples_medians.tsv -p MUSMU -o MUSMU_tidy_corr.csv 
```
Oyster
```
python3 /proj/robert/scripts/density/block_correlation_analysis.py -t 02_gene_density_analysis/density_whole_genome/key_nodes.tidydf.csv -s 01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -e CRAGI/CRAGI_transcript_tpms_all_samples.tsv -p CRAGI -o CRAGI_tidy_corr.csv 
```
Urchin
```
python3 /proj/robert/scripts/density/block_correlation_analysis.py -t 02_gene_density_analysis/density_whole_genome/key_nodes.tidydf.csv -s 01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -e STRPU/STRPU_transcript_tpms_all_samples.tsv -p STRPU -o STRPU_tidy_corr.csv 
```
Hemichordate
```
python3 /proj/robert/scripts/density/block_correlation_analysis.py -t 02_gene_density_analysis/density_whole_genome/key_nodes.tidydf.csv -s 01_microsynteny/chrom_of/5.blocks.3.syn_corrected.synt -e SACKO/SACKO_transcript_tpms_all_samples.tsv -p SACKO -o SACKO_tidy_corr.csv 
```
## 8. Particular cases manual annotation
### 8.1 HOX
#### 8.1.1 HOX annotation
References and animals for which *hox* sequences used as queries for reciprocal BLAST searches in the planulozoan sequences of our sample:

* Amemiya et al. 2013 (doi: doi.org/10.1038/nature12027); LATCH

* Anaya et al. 2013 (doi:10.1186/1471-213X-13-26); HOMSA, MUSMU, DANRE, CALMI, CHEMY, CIOIN, BRALA, STRPU, SACKO, PTYFL

* Belcaid et al. 2018 (doi: 10.1073/pnas.1817322116); EUPSC

* Brauchle et al. 2018 (doi: doi.org/10.1093/gbe/evy170); HOFMI

* Currie et al. 2016 (doi: 10.1186/s13227-016-0044-8); SCHME

* DuBuc et al. 2012 (doi: 10.1093/icb/ics098); NEMVE

* Hench et al. 2015 (doi: 10.1371/journal.pone.0126947); CAEEL

* Kuraku and Meyer 2009 (doi: 10.1387/ijdb.072533km); DANRE, HOMSA, GALGA, XENTR, CALMI

* Leclère et al. 2019 (doi: 10.1038/s41559-019-0833-2); CLYHE

* Leite et al. 2018 (doi: 10.1093/molbev/msy125); PARTE

* Pace et al. 2016 (doi: 10.1186/s13227-016-0048-4); ANOGA, DROME, TRICA, DAPPU, STRMA, IXOSC

* Simakov et al. 2013 (doi: 10.1038/nature11696); CAPTE, LOTGI, DAPPU

* Simakov et al. 2015 (doi: 10.1038/nature16150); PTYFL, SACKO

* Wang et al. 2017 (doi: 10.1038/s41559-017-0120); MIZYE

Orthologs of hoxes in the 49 species of our sample were identified by reciprocal BLASTs using as queries aforementionned sequences.
Annotation was given according to their top hit in the NR database.

#### 8.1.2 Hox microsynteny and gene density analysis
`HOX.clus` Is a cluster file constituted of only one OG (accessions of all annoatted hoxes in http://synteny.csb.univie.ac.at/, planulozoan novelty, `hox manually curated`).

```
cd 01_microsynteny/chrom
prepMicroSynt.pl ACAPL.chrom,ACRMI.chrom,ADIVA.chrom,AMPQU.chrom,ANOGA.chrom,AURAU.chrom,BRALA.chrom,CAEEL.chrom,CALMI.chrom,CAPOW.chrom,CAPTE.chrom,CHEMY.chrom,CIOIN.chrom,CLYHE.chrom,CRAGI.chrom,DANRE.chrom,DAPPU.chrom,DROME.chrom,EUPSC.chrom,EXAPA.chrom,GALGA.chrom,HELRO.chrom,HIPCO.chrom,HOFMI.chrom,HOIHO.chrom,HOMSA.chrom,HYDVU.chrom,IXOSC.chrom,LATCH.chrom,LEPOC.chrom,LINAN.chrom,LOTGI.chrom,MAYZE.chrom,MIZYE.chrom,MNELE.chrom,MUSMU.chrom,NEMVE.chrom,PARTE.chrom,PLEBA.chrom,PTYFL.chrom,SACKO.chrom,SALRO.chrom,SCHME.chrom,STRMA.chrom,STRPU.chrom,SYCCI.chrom,TRIAD.chrom,TRICA.chrom,XENTR.chrom 5 03_HOX/HOX_density/HOX.clus
sbatch --array=1-1176 --constraint=array-1core --job-name=hox_synt job.sh

chromfile_folder='01_microsynteny/chrom'
makeClusters3.pl $chromfile_folder/ACAPL.chrom,$chromfile_folder/ACRMI.chrom,$chromfile_folder/ADIVA.chrom,$chromfile_folder/AMPQU.chrom,$chromfile_folder/ANOGA.chrom,$chromfile_folder/AURAU.chrom,$chromfile_folder/BRALA.chrom,$chromfile_folder/CAEEL.chrom,$chromfile_folder/CALMI.chrom,$chromfile_folder/CAPOW.chrom,$chromfile_folder/CAPTE.chrom,$chromfile_folder/CHEMY.chrom,$chromfile_folder/CIOIN.chrom,$chromfile_folder/CLYHE.chrom,$chromfile_folder/CRAGI.chrom,$chromfile_folder/DANRE.chrom,$chromfile_folder/DAPPU.chrom,$chromfile_folder/DROME.chrom,$chromfile_folder/EUPSC.chrom,$chromfile_folder/EXAPA.chrom,$chromfile_folder/GALGA.chrom,$chromfile_folder/HELRO.chrom,$chromfile_folder/HIPCO.chrom,$chromfile_folder/HOFMI.chrom,$chromfile_folder/HOIHO.chrom,$chromfile_folder/HOMSA.chrom,$chromfile_folder/HYDVU.chrom,$chromfile_folder/IXOSC.chrom,$chromfile_folder/LATCH.chrom,$chromfile_folder/LEPOC.chrom,$chromfile_folder/LINAN.chrom,$chromfile_folder/LOTGI.chrom,$chromfile_folder/MAYZE.chrom,$chromfile_folder/MIZYE.chrom,$chromfile_folder/MNELE.chrom,$chromfile_folder/MUSMU.chrom,$chromfile_folder/NEMVE.chrom,$chromfile_folder/PARTE.chrom,$chromfile_folder/PLEBA.chrom,$chromfile_folder/PTYFL.chrom,$chromfile_folder/SACKO.chrom,$chromfile_folder/SALRO.chrom,$chromfile_folder/SCHME.chrom,$chromfile_folder/STRMA.chrom,$chromfile_folder/STRPU.chrom,$chromfile_folder/SYCCI.chrom,$chromfile_folder/TRIAD.chrom,$chromfile_folder/TRICA.chrom,$chromfile_folder/XENTR.chrom .5.blocks 3 0.3 0.5  03_HOX/HOX_density/5.blocks.3.syn.synt

correct_blocks_coordinates.py 03_HOX/HOX_density/5.blocks.3.syn.synt ACAPL.chrom,ACRMI.chrom,ADIVA.chrom,AMPQU.chrom,ANOGA.chrom,AURAU.chrom,BRALA.chrom,CAEEL.chrom,CALMI.chrom,CAPOW.chrom,CAPTE.chrom,CHEMY.chrom,CIOIN.chrom,CLYHE.chrom,CRAGI.chrom,DANRE.chrom,DAPPU.chrom,DROME.chrom,EUPSC.chrom,EXAPA.chrom,GALGA.chrom,HELRO.chrom,HIPCO.chrom,HOFMI.chrom,HOIHO.chrom,HOMSA.chrom,HYDVU.chrom,IXOSC.chrom,LATCH.chrom,LEPOC.chrom,LINAN.chrom,LOTGI.chrom,MAYZE.chrom,MIZYE.chrom,MNELE.chrom,MUSMU.chrom,NEMVE.chrom,PARTE.chrom,PLEBA.chrom,PTYFL.chrom,SACKO.chrom,SALRO.chrom,SCHME.chrom,STRMA.chrom,STRPU.chrom,SYCCI.chrom,TRIAD.chrom,TRICA.chrom,XENTR.chrom > 03_HOX/HOX_density/5.blocks.3.syn_corrected.synt

```
Make the randomized HOX blocks.
```
cd 03_HOX/HOX_density/rand_blocks
cut -f2 ../5.blocks.3.syn_corrected.synt|sort -u > genomes.list
n=100
outfolder='03_HOX/HOX_density/rand_blocks'
chromfile_folder='01_microsynteny/chrom'
inputfile='../5.blocks.3.syn_corrected.synt'

while read species; do echo "pick_random_blocks.R --only-genome=$species -n $n $outfolder $inputfile $chromfile_folder/$species.chrom"; done < genomes.list > random_hoxblocks

slurmtasks -m 10 -n random_hoxblocks random_hoxblocks | sbatch
cat *.tsv | grep -v block_id > Hoxblocks.100r.syntx #merge all the tsvs together

cd $chromfile_folder
correct_blocks_coordinates.py 03_HOX/HOX_density/rand_blocks/Hoxblocks.100r.syntx ACAPL.chrom,ACRMI.chrom,ADIVA.chrom,AMPQU.chrom,ANOGA.chrom,AURAU.chrom,BRALA.chrom,CAEEL.chrom,CALMI.chrom,CAPOW.chrom,CAPTE.chrom,CHEMY.chrom,CIOIN.chrom,CLYHE.chrom,CRAGI.chrom,DANRE.chrom,DAPPU.chrom,DROME.chrom,EUPSC.chrom,EXAPA.chrom,GALGA.chrom,HELRO.chrom,HIPCO.chrom,HOFMI.chrom,HOIHO.chrom,HOMSA.chrom,HYDVU.chrom,IXOSC.chrom,LATCH.chrom,LEPOC.chrom,LINAN.chrom,LOTGI.chrom,MAYZE.chrom,MIZYE.chrom,MNELE.chrom,MUSMU.chrom,NEMVE.chrom,PARTE.chrom,PLEBA.chrom,PTYFL.chrom,SACKO.chrom,SALRO.chrom,SCHME.chrom,STRMA.chrom,STRPU.chrom,SYCCI.chrom,TRIAD.chrom,TRICA.chrom,XENTR.chrom > 03_HOX/HOX_density/rand_blocks/Hoxblocks.100r.synt

```

```
python3 make_tidy_density_df.py -s ../../03_REDUX_density_HOX_WNT/hox -g ../../02_REDUX_gene_density_analysis/genomes/ -c ../../01_microsynteny/chrom/ -m ../../03_HOX/hox_graph_data/hox_planu.clus -og ../..//03_HOX/hox_graph_data/hox_planu.clus -o hox
```


### 8.2 WNT
#### 8.2.1 Characterization of wnt5-wnt7 cluster
The cluster detected by the microsynteny pipeline is cluster 953, planulozoan novelty (http://synteny.csb.univie.ac.at/).

References for evidence of synteny of members of the wnt5-wnt7 cluster:

* Cho et al. 2011 (doi: 10.1093/molbev/msq052)

* Janssen et al. 2010 (doi: 10.1186/1471-2148-10-374)

* Sullivan et al. 2007 (doi: 10.1007/s00427-007-0136-5)

* Irimia et al. 2012 (doi: 10.1101/gr.139725.112)

* Kapasa et al. 2010 (doi: 10.1186/1745-6150-5-49)

* Garriock et al. 2012 (doi: 10.1002/dvdy.21156)

Known syntenies by reference

* *wnt5* and *wnt7* in Planulozoans: (Cho et al. 2011), (Janssen et al. 2010), (Sullivan et al. 2007)

* *wnt5* and *fbxl14* in Planulozoans (Irimia et al. 2012), and bilaterians (Kapasa et al. 2010)

* *wnt7* and *atxn10* in Planulozoans (Irimia et al. 2012) in Vertebrates (Garriock et al. 2012)

* *wnt5* and *erc1/2* in Olfactores (Kapasa at al. 2010)

* *wnt5* and *cacna2d* in Vertebrates (Kapasa et al. 2010, Garriock et al. 2012)

* *cacna1d*, *ninj1/2* and *dcp1*, this study


Orthologs to members of the *wnt5-wnt7* cluster (*wnt5-wnt7* pair, *fbxl14*, *atxn10*, *erc1/2*, *cacna1d*, *cacna2d* *ninj1/ninj2* and *dcp1a/dcp1b*) in the planulozoans of the 49 species of our sample were identified by reciprocal BLAST using as inital queries human sequences.
#### 8.2.2 wnt5-wnt7 microsynteny and gene density analysis

Similarmy to what we did with hoxes, we build an orthology file where all members of the cluster are in a single OG (*wnt5-wnt7* pair, *fbxl14*, *atxn10*, *erc1/2*, *cacna1d*, *cacna2d* *ninj1/ninj2* and *dcp1a/dcp1b*), the file is `wnt5_wnt7.clus`

```
cd 01_microsynteny/chrom
perl /proj/Simakov/scripts/MICROSYNT/prepMicroSynt.pl ACAPL.chrom,ACRMI.chrom,ADIVA.chrom,ANOGA.chrom,AURAU.chrom,BRALA.chrom,CAEEL.chrom,CALMI.chrom,CAPTE.chrom,CHEMY.chrom,CIOIN.chrom,CLYHE.chrom,CRAGI.chrom,DANRE.chrom,DAPPU.chrom,DROME.chrom,EUPSC.chrom,EXAPA.chrom,GALGA.chrom,HELRO.chrom,HIPCO.chrom,HOFMI.chrom,HOMSA.chrom,HYDVU.chrom,IXOSC.chrom,LATCH.chrom,LEPOC.chrom,LINAN.chrom,LOTGI.chrom,MAYZE.chrom,MIZYE.chrom,MUSMU.chrom,NEMVE.chrom,PARTE.chrom,PTYFL.chrom,SACKO.chrom,SCHME.chrom,STRMA.chrom,STRPU.chrom,TRICA.chrom,XENTR.chrom 5 07_wnt5_wnt7/wnt5_wnt7.clus
sbatch --array=1-820 --constraint=array-1core --job-name=wnt5_wnt7 job.sh

chromfile_folder='01_microsynteny/chrom'
perl /proj/Simakov/scripts/MICROSYNT/makeClusters3.pl $chromfile_folder/ACAPL.chrom,$chromfile_folder/ACRMI.chrom,$chromfile_folder/ADIVA.chrom,$chromfile_folder/AMPQU.chrom,$chromfile_folder/ANOGA.chrom,$chromfile_folder/AURAU.chrom,$chromfile_folder/BRALA.chrom,$chromfile_folder/CAEEL.chrom,$chromfile_folder/CALMI.chrom,$chromfile_folder/CAPOW.chrom,$chromfile_folder/CAPTE.chrom,$chromfile_folder/CHEMY.chrom,$chromfile_folder/CIOIN.chrom,$chromfile_folder/CLYHE.chrom,$chromfile_folder/CRAGI.chrom,$chromfile_folder/DANRE.chrom,$chromfile_folder/DAPPU.chrom,$chromfile_folder/DROME.chrom,$chromfile_folder/EUPSC.chrom,$chromfile_folder/EXAPA.chrom,$chromfile_folder/GALGA.chrom,$chromfile_folder/HELRO.chrom,$chromfile_folder/HIPCO.chrom,$chromfile_folder/HOFMI.chrom,$chromfile_folder/HOIHO.chrom,$chromfile_folder/HOMSA.chrom,$chromfile_folder/HYDVU.chrom,$chromfile_folder/IXOSC.chrom,$chromfile_folder/LATCH.chrom,$chromfile_folder/LEPOC.chrom,$chromfile_folder/LINAN.chrom,$chromfile_folder/LOTGI.chrom,$chromfile_folder/MAYZE.chrom,$chromfile_folder/MIZYE.chrom,$chromfile_folder/MNELE.chrom,$chromfile_folder/MUSMU.chrom,$chromfile_folder/NEMVE.chrom,$chromfile_folder/PARTE.chrom,$chromfile_folder/PLEBA.chrom,$chromfile_folder/PTYFL.chrom,$chromfile_folder/SACKO.chrom,$chromfile_folder/SALRO.chrom,$chromfile_folder/SCHME.chrom,$chromfile_folder/STRMA.chrom,$chromfile_folder/STRPU.chrom,$chromfile_folder/SYCCI.chrom,$chromfile_folder/TRIAD.chrom,$chromfile_folder/TRICA.chrom,$chromfile_folder/XENTR.chrom .5.blocks 3 0.3 0.5 > 07_wnt5_wnt7/wnt5_wnt7_density/pairwise_blocks/5.blocks.3.syn.synt

cd $chromfile_folder
correct_blocks_coordinates.py 07_wnt5_wnt7/wnt5_wnt7_density/5.blocks.3.syn.synt ACAPL.chrom,ACRMI.chrom,ADIVA.chrom,AMPQU.chrom,ANOGA.chrom,AURAU.chrom,BRALA.chrom,CAEEL.chrom,CALMI.chrom,CAPOW.chrom,CAPTE.chrom,CHEMY.chrom,CIOIN.chrom,CLYHE.chrom,CRAGI.chrom,DANRE.chrom,DAPPU.chrom,DROME.chrom,EUPSC.chrom,EXAPA.chrom,GALGA.chrom,HELRO.chrom,HIPCO.chrom,HOFMI.chrom,HOIHO.chrom,HOMSA.chrom,HYDVU.chrom,IXOSC.chrom,LATCH.chrom,LEPOC.chrom,LINAN.chrom,LOTGI.chrom,MAYZE.chrom,MIZYE.chrom,MNELE.chrom,MUSMU.chrom,NEMVE.chrom,PARTE.chrom,PLEBA.chrom,PTYFL.chrom,SACKO.chrom,SALRO.chrom,SCHME.chrom,STRMA.chrom,STRPU.chrom,SYCCI.chrom,TRIAD.chrom,TRICA.chrom,XENTR.chrom > 07_wnt5_wnt7/wnt5_wnt7_density/5.blocks.3.syn_corrected.synt

cd 07_wnt5_wnt7/wnt5_wnt7_density/rand_blocks

cut -f2 ../5.blocks.3.syn_corrected.synt|sort -u > genomes.list
n=100
outfolder='07_wnt5_wnt7/wnt5_wnt7_density/rand_blocks'
chromfile_folder='01_microsynteny/chrom'
inputfile='../5.blocks.3.syn_corrected.synt'

while read species; do echo "Rscript  /proj/robert/scripts/microsynteny/pick_random_blocks.R --only-genome=$species -n $n $outfolder $inputfile $chromfile_folder/$species.chrom"; done < genomes.list > wnt5_wnt7

nohup bash wnt5_wnt7&

cat *.tsv | grep -v block_id > wnt5_wnt7.100r.syntx #merge all the tsvs together

cd $chromfile_folder
correct_blocks_coordinates.py 07_wnt5_wnt7/wnt5_wnt7_density/rand_blocks/wnt5_wnt7.100r.syntx ACAPL.chrom,ACRMI.chrom,ADIVA.chrom,AMPQU.chrom,ANOGA.chrom,AURAU.chrom,BRALA.chrom,CAEEL.chrom,CALMI.chrom,CAPOW.chrom,CAPTE.chrom,CHEMY.chrom,CIOIN.chrom,CLYHE.chrom,CRAGI.chrom,DANRE.chrom,DAPPU.chrom,DROME.chrom,EUPSC.chrom,EXAPA.chrom,GALGA.chrom,HELRO.chrom,HIPCO.chrom,HOFMI.chrom,HOIHO.chrom,HOMSA.chrom,HYDVU.chrom,IXOSC.chrom,LATCH.chrom,LEPOC.chrom,LINAN.chrom,LOTGI.chrom,MAYZE.chrom,MIZYE.chrom,MNELE.chrom,MUSMU.chrom,NEMVE.chrom,PARTE.chrom,PLEBA.chrom,PTYFL.chrom,SACKO.chrom,SALRO.chrom,SCHME.chrom,STRMA.chrom,STRPU.chrom,SYCCI.chrom,TRIAD.chrom,TRICA.chrom,XENTR.chrom > 07_wnt5_wnt7/wnt5_wnt7_density/rand_blocks/wnt5_wnt7.100r.synt


python3 make_tidy_density_df.py -s ../../03_HOX_WNT/wnt -g ../../02_REDUX_gene_density_analysis/genomes/ -c ../../01_microsynteny/chrom/ -m ../../07_wnt5_wnt7/NeuroWnt_density_manual_redone/manual_final.clusters -og /scratch/robert/2019_microsynteny_size_constraints/07_wnt5_wnt7/NeuroWnt_density_manual_redone/wnt5_7_graph_data/wnt5_wnt7.clus -o wnt

```


Regarding blocks, some of them a lineage specific expansions of CACNA1 or CACNA2D subunits. We filter the blocks for keeping only the ones with *wnt5* or *wnt7*## 9 Graph representation of blocks
## 2. File preparation
get assembly length of the genomes. In the normali (different from lengths than the used for normalizing)
```
for i in $(ls 02_gene_density_analysis/genomes/*genome); do get_length_genome.py -n 0; done
```
Convert the blocks to edgelists (graphs where OGs are the nodes, and edges are the distances of syntenic orthogroups.
We also make edgelists for manually curated hox and wnt genes.
```
Block_to_OGcommus.py -c 02_gene_density_analysis/gene_density/second_sampling/Metazoa_total/Metazoa.m2.total.clusters -b 02_gene_density_analysis/gene_density/second_sampling/Metazoa_total/Metazoa.m2.total.synt -g ACAPL.chrom ACRMI.chrom ADIVA.chrom AMPQU.chrom ANOGA.chrom AURAU.chrom BRALA.chrom CAEEL.chrom CALMI.chrom CAPOW.chrom CAPTE.chrom CIOIN.chrom CLYHE.chrom CRAGI.chrom DAPPU.chrom DROME.chrom EUPSC.chrom EXAPA.chrom HELRO.chrom HOFMI.chrom HOIHO.chrom HYDVU.chrom IXOSC.chrom LATCH.chrom LINAN.chrom LOTGI.chrom MIZYE.chrom MNELE.chrom NEMVE.chrom PARTE.chrom PLEBA.chrom PTYFL.chrom SACKO.chrom SCHME.chrom STRMA.chrom STRPU.chrom SYCCI.chrom TRIAD.chrom TRICA.chrom XENTR.chrom -og 01_microsynteny/Orthofinder.clus -o 10_OG_graphs_representation/Metazoa_total -r 02_gene_density_analysis/gene_density/second_sampling/Metazoa_total/rand_blocks/Metazoa.m2.total.100r.synt

Block_to_OGcommus.py -c 02_gene_density_analysis/gene_density/second_sampling/Planu_placo/Planu_Placo.m2.novel.clusters -b 02_gene_density_analysis/gene_density/second_sampling/Planu_placo/Planu_Placo.m2.novel.synt -g ACRMI.chrom ADIVA.chrom ANOGA.chrom BRALA.chrom CAPTE.chrom CIOIN.chrom CRAGI.chrom EUPSC.chrom EXAPA.chrom HOIHO.chrom HYDVU.chrom LINAN.chrom LOTGI.chrom MIZYE.chrom PARTE.chrom PTYFL.chrom SACKO.chrom SCHME.chrom STRPU.chrom TRIAD.chrom TRICA.chrom -og 01_microsynteny/Orthofinder.clus -o 10_OG_graphs_representation/Planu_Placo_novel -r 02_gene_density_analysis/gene_density/second_sampling/Planu_placo/rand_blocks/Planu_Placo.m2.novel.100r.synt

Block_to_OGcommus.py -c 02_gene_density_analysis/gene_density/second_sampling/Bilateria/Bilateria.m2.novel.clusters -b 02_gene_density_analysis/gene_density/second_sampling/Bilateria/Bilateria.m2.novel.synt -g ACAPL.chrom ADIVA.chrom ANOGA.chrom BRALA.chrom CAEEL.chrom CALMI.chrom CAPTE.chrom CHEMY.chrom CIOIN.chrom CRAGI.chrom DANRE.chrom DAPPU.chrom DROME.chrom EUPSC.chrom GALGA.chrom HELRO.chrom HIPCO.chrom HOFMI.chrom HOMSA.chrom IXOSC.chrom LATCH.chrom LEPOC.chrom LINAN.chrom LOTGI.chrom MAYZE.chrom MIZYE.chrom MUSMU.chrom PARTE.chrom PTYFL.chrom SACKO.chrom SCHME.chrom STRMA.chrom STRPU.chrom TRICA.chrom XENTR.chrom -og 01_microsynteny/Orthofinder.clus -o 10_OG_graphs_representation/Bilateria_novel -r 02_gene_density_analysis/gene_density/second_sampling/Bilateria/rand_blocks/Bilateria.m2.novel.100r.synt

Block_to_OGcommus.py -c 02_gene_density_analysis/gene_density/second_sampling/Planulozoa/Planulozoa.m2.novel.clusters -b 02_gene_density_analysis/gene_density/second_sampling/Planulozoa/Planulozoa.m2.novel.synt -g ACAPL.chrom ACRMI.chrom ADIVA.chrom ANOGA.chrom AURAU.chrom BRALA.chrom CAEEL.chrom CALMI.chrom CAPTE.chrom CHEMY.chrom CIOIN.chrom CLYHE.chrom CRAGI.chrom DANRE.chrom DAPPU.chrom DROME.chrom EUPSC.chrom EXAPA.chrom GALGA.chrom HELRO.chrom HIPCO.chrom HOFMI.chrom HOMSA.chrom HYDVU.chrom IXOSC.chrom LATCH.chrom LEPOC.chrom LINAN.chrom LOTGI.chrom MAYZE.chrom MIZYE.chrom MUSMU.chrom NEMVE.chrom PARTE.chrom PTYFL.chrom SACKO.chrom SCHME.chrom STRMA.chrom STRPU.chrom TRICA.chrom XENTR.chrom -og 01_microsynteny/Orthofinder.clus -o 10_OG_graphs_representation/Planulozoa_novel -r 02_gene_density_analysis/gene_density/second_sampling/Planulozoa/rand_blocks/Planulozoa.m2.novel.100r.synt

Block_to_OGcommus.py -c 02_gene_density_analysis/gene_density/second_sampling/recent_blocks/Lophotrochozoa/Lophotrochozoa_novel.clusters -b 02_gene_density_analysis/gene_density/second_sampling/recent_blocks/Lophotrochozoa/Lophotrochozoa_novel.synt -g ADIVA.chrom CAPTE.chrom CRAGI.chrom EUPSC.chrom HELRO.chrom LINAN.chrom LOTGI.chrom MIZYE.chrom SCHME.chrom -og 01_microsynteny/Orthofinder.clus -o 10_OG_graphs_representation/Lophotrochozoa_novel -r 02_gene_density_analysis/gene_density/second_sampling/recent_blocks/Lophotrochozoa/rand_blocks/Lophotrochozoa_novel.100r.synt

Block_to_OGcommus.py -c 02_gene_density_analysis/gene_density/second_sampling/recent_blocks/Vertebrata/Vertebrata_novel.clusters -b 02_gene_density_analysis/gene_density/second_sampling/recent_blocks/Vertebrata/Vertebrata_novel.synt -g CALMI.chrom CHEMY.chrom DANRE.chrom GALGA.chrom HIPCO.chrom HOMSA.chrom LATCH.chrom LEPOC.chrom MAYZE.chrom MUSMU.chrom XENTR.chrom -og 01_microsynteny/Orthofinder.clus -o 10_OG_graphs_representation/Vertebrata_novel -r 02_gene_density_analysis/gene_density/second_sampling/recent_blocks/Vertebrata/rand_blocks/Vertebrata_novel.100r.synt


```

We use the lengths obtained with `get_length_genome.py` to normalise the basepairs distances in the OG pairs edgelists.


WNT graph
We'll use our annotation. There is 8 OGs that are annotated as being part of the syntenic block. We'll have 8 OGs in the 07_wnt5_wnt7/wnt5_7_graph_data/wnt5_wnt7.clus.

ATXN10, CACNA1, CACNA2,DCP1, ERC, FBXL14, NINJ, WNT
```
cd 01_microsynteny/chrom
Block_to_OGcommus_2.py -c 07_wnt5_wnt7/manual_final.clusters -b 07_wnt5_wnt7/manual_final.synt -g ACAPL.chrom ACRMI.chrom ANOGA.chrom BRALA.chrom CALMI.chrom CHEMY.chrom DANRE.chrom DAPPU.chrom EUPSC.chrom EXAPA.chrom GALGA.chrom HIPCO.chrom HOMSA.chrom IXOSC.chrom LATCH.chrom LEPOC.chrom LINAN.chrom LOTGI.chrom MAYZE.chrom MUSMU.chrom NEMVE.chrom PARTE.chrom SACKO.chrom TRICA.chrom XENTR.chrom  -og 07_wnt5_wnt7/wnt5_7_graph_data/wnt5_wnt7.clus -o 07_wnt5_wnt7/wnt5_7_graph_data/wnt5_wnt7 -r 07_wnt5_wnt7/rand_blocks/manual_final.100r.synt --custom_orthology
```

HOX graph
We'll use our annotation `03_HOX/hox_graph_data/hox_planu.clus`, where hoxes are grouped by taxonomic groups (Lophotrochozoan hox OGs, Vertebrate hox Ogs, etc.) specific.
```
cd 01_microsynteny/chrom
Block_to_OGcommus_2.py -c 03_HOX/HOX_density/5.blocks.3.syn.clusters -b 03_HOX/HOX_density/5.blocks.3.syn_corrected.synt -g ACAPL.chrom ACRMI.chrom ANOGA.chrom BRALA.chrom CALMI.chrom CAPTE.chrom CHEMY.chrom CIOIN.chrom CRAGI.chrom DAPPU.chrom DROME.chrom EUPSC.chrom EXAPA.chrom GALGA.chrom HELRO.chrom HIPCO.chrom HOMSA.chrom IXOSC.chrom LATCH.chrom LEPOC.chrom LINAN.chrom LOTGI.chrom MAYZE.chrom MIZYE.chrom MUSMU.chrom NEMVE.chrom PARTE.chrom PTYFL.chrom SACKO.chrom STRMA.chrom STRPU.chrom TRICA.chrom XENTR.chrom -og 03_HOX/hox_graph_data/hox_planu.clus -o 03_HOX/hox_graph_data/hox -r 03_HOX/HOX_density/rand_blocks/Hoxblocks.100r.synt --custom_orthology
```

