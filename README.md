# **LOCATOR: A toolkit for incorporating local ancestry in association tests**

LOCATOR is a comprehensive toolkit for analyzing and incorporating local ancestry information in Genome-wide Association Studies (GWAS). This package provides tools for manipulating local ancestry estimates and accounting for local ancestry effects in association testing.
This tutorial covers:
  - Output data structures
  - Data preparation guidelines
  - Pipeline construction
  - Caveats

For detailed technical specifications and methodology, please refer to LOCATOR package manual.pdf and *A Continuous Local Ancestry Measure for Efficient Local-Ancestry-Aware Association Tests. Hanxiao Sun. (2024)* https://digitalcommons.library.tmc.edu/uthsph_dissertsopen/248/.

Getting Help<br>
If you have questions or encounter issues, please open a discussion thread in this repository. We actively monitor and respond to user inquiries.

Citation<br>
If you use LOCATOR in your research, please cite: Sun, Hanxiao, "A Continuous Local Ancestry Measure for Efficient Local-Ancestry-Aware Association Tests" (2024). Dissertations & Theses (Open Access). 248. https://digitalcommons.library.tmc.edu/uthsph_dissertsopen/248.

I would also like to credit this project to Dr. Han Chen, Shuyi Guo from University of Texas Health Science Center, Tomin Perea-Chamblee from Memorial Sloan Kettering Cancer Center and other authors in developing ideas, finalizing results and preparing the whole R package.
A detailed instruction and manuscript will be available soon.

## Rationale
Local Ancestry Coordinates (LACs) are derived from the established linear relationships between genetic ancestry estimates and Principal Components (PCs). We first establish associations between global ancestry and global PCs to identify projection anchors representing hypothetical homogeneous ancestry populations. Since local ancestry inferences can be interpreted as ancestry composition probabilities, this relationship allows us to project these inferences onto the space spanned by global PCs, with the projections weighted by local ancestry probabilities.
Let $`X_g`$, $`X_l`$ as global PCs and LACs, $`W_g`$, $`W_l`$ as global and local ancestry, $`\psi`$ as the ancestry projection matrix. There is

$$X_g = W_g\psi+\epsilon$$
$$\hat{\psi} = (W_g^TW_g)^{-1}W_g^TX_g$$
$$\hat{X_l} = W_l(W_g^TW_g)^{-1}W_g^TX_g$$

LACs serve as efficient proxies for local ancestry estimates in GWAS adjustment. Due to memory constraints, these coordinates are stored in binary format. For the complete mathematical framework underlying this approach, please refer to *A Continuous Local Ancestry Measure for Efficient Local-Ancestry-Aware Association Tests. Hanxiao Sun. (2024)*. Detailed information about the data structure can be found in the section Step 2.5.

## Example data
We simulated a group of 200 three-way admixed Hispanic individuals from African (AFR), European (EUR), and Amerindian (AMR) by Admix-Simu and selected a 20 Mbp region on chromosome 12 (chr12:48275834 - 68009712) containing 45,170 SNPs with minor allele frequency (MAF) ≥ 0.05 and 100 local ancestry tracks. Please note that this example dataset was simulated for demonstration purposes only, and the local ancestry was inferred using a combined continental reference panel, which may introduce some imprecision and does not fully capture the genetic complexity of real admixed populations.

### Optional: Reference panel preparation
A reference panel that includes the ancestral components present in your test samples is recommended, though not strictly required. If needed, please prepare your reference panel to match your study population, plus additional relevant reference populations to ensure accurate ancestral representation and positioning if applicable. To establish a consistent ancestry representation framework, please first perform Principal Component Analysis (PCA) on the reference panel, then project your test samples onto this reference-panel-defined space. Alternatively, please provide the self-adjusted global PCs instead.

### Step 0: Data preparation, extract local ancestry breakpoint and calculate global ancestry
As genetic variants within the same ancestry tract share identical local ancestry assignments, we first use get.breakpoints() to identify and extract local ancestry breakpoints from the input local ancestry inference. Before using get.breakpoints(), please format the input local ancestry data properly. The toolkit accepts discrete local ancestry inferences from any software platform. For RFMix v2.0 users, please remove the first header row and separate the first four columns (genetic position information) from the local ancestry inference section. Please ensure local ancestry data follows this order: hap1:::ancestry1, hap1:::ancestry2, hap2:::ancestry1, hap2:::ancestry2, with consistent ancestry labels throughout.

The get.breakpoints() function simultaneously identifies local ancestry breakpoints and calculates local-ancestry-aggregated global ancestry estimates. Please note that get.breakpoints() utilizes Linux/Unix commands to read large input files in chunks (1000 rows per time), so please grant R appropriate file access permissions.

We also provide internal local-ancestry-aggregated global ancestry estimates with genetic positions counted from telomeres (start.pos = 0). Please toggle on the argument if.yield.ga = TRUE. For analysis targeted at specific SNP chunks, please adjust start.pos accordingly. Alternatively, external global ancestry estimates (e.g., from ADMIXTURE) are also accepted for downstream analysis. Please ensure it is consistent with local ancestry inference in this case.

```
### Generate breakpoint information files
### Please be aware that we currently don't provide example dataset for local ancestry inference in terms of the file size. This chunk of code is used for demonstration only.
get.breakpoints(input.la=%placeholder%,input.pos=%placeholder%,input.id=%placeholder%,input.anc=c('AFR','EUR','AMR'),input.sep='\t',if.yield.ga=T,output.bp=%placeholder%,output.ga=%placeholder%)  ## please replace all placeholders with your file paths and parameter settings.
```

### Step 1: Calculate the ancestry projection matrix $`\psi`$
The ancestry projection matrix $`\psi`$ is calculated using the get.anchor() function, which provides the exact positioning information for hypothetically unadmixed populations. Please only include the dimensions of global PCs that will be utilized in the downstream analysis as needed. Please ensure the order of global ancestry matches the local ancestry inference that will be analyzed in subsequent steps.

```
ga_file_path<-system.file('estdata','global_ancestry',package='LOCATOR')
ga_file<-fread(ga_file_path)
gc_file_path<-system.file('estdata','geno_global_PCs',package='LOCATOR')
gc_file<-fread(gc_file_path)
psi<-get.anchor(input.ga=ga_file[,2:4],input.gc=gc_file[,2:6])
```

### Step 2: Calculate Local Ancestry Coordinates (LACs)
Here we provide two types of LACs transformed from extracted local ancestry breakpoints: original LACs and orthogonalized LACs, the latter of which are derived by removing global PC effects from original LACs through orthogonal decomposition. We define the difference between orthogonalized LACs and global PCs as Local Deviation $`D_L`$, which represents pure additional population stratification effects derived from local ancestry. Please choose the LAC type that best fits your analysis needs. Please generate original LACs via get.LACs() in advance as the required input for orthogonalized LACs via get.orthogonal.LACs(). When applying these functions, please specify the type of sample IDs, selected IDs, global PCs, positions as needed. Both LAC types are stored in the same binary format as breakpoints. Let $`L`$ be the orthogonalized LACs, we have:

$$L = \[I_n - X_g(X_g^TX_g)^{-1}X_g\]W_l\hat{\psi}$$
$$D_L = X_g - L $$

```
### Generate original LACs files
breakpoint_file_path<-system.file('estdata','breakpoint_info',package='LOCATOR')
get.LACs(input.bp=breakpoint_file_path,psi=psi,output.la=%placeholder%,id.type='character')  ## please put down your LACs output filename here
### Generate orthogonalized LACs files
### Please first convert sample IDs into rownames
gc_file<-gc_file %>% tibble::column_to_rownames('IID')
get.orthogonal.LACs(input.la='tmp5_LACs',input.gc=gc_file,output.la=%placeholder%,id.type='character')  ## please put down your LACs output filename here
```

### Step 2.5 (Optional): View binary files
Please be aware that local ancestry breakpoints and LACs are stored in binary format and should not be opened or modified directly to prevent potential file corruption. Instead, please use the read.binary() function to view these binary files. If you encounter any file damage (e.g.: missing values), please examine your input data and regenerate the binary files.
The binary files constitute two sections. The first section is about file information. In each file, it is arranged in the order of:
  - The number of PCs involved
  - The number of samples documented
  - The list of exact sample IDs
  - The number of genetic positions documented
  - The list of exact genetic positions
where sample IDs and genetic positions are stored as vectors.
The second section is about actual breakpoints or LAC values. Before calling read.binary(), please specify the file type, the sample ID type, and selected positions and PCs. It will return a list of queries.

```
### View local ancestry information at the first breakpoint position
bp_slot<-read.binary(input.file=breakpoint_file_path,type='breakpoints',id.type='character',select.pos=48275834)
View(bp_slot[[1]])

sampleIDs	AFR	EUR	AMR
Admixed2	1.0	0.0	0.0
Admixed3	0.5	0.0	0.5
Admixed7	1.0	0.0	0.0
Admixed9	0.0	0.0	1.0
Admixed11	1.0	0.0	0.0
Admixed18	0.5	0.0	0.5
```

### Step 3: Perform local-ancestry-aware association test
We use get.la.aware.score() to perform local-ancestry-aware association tests through a two-step process: First, please construct a null model using glmmkin() from the GMMAT R package to adjust for covariates, including both fixed and random effects. Please refer to [https://github.com/hanchenphd/GMMAT] for detailed instructions on null model construction. Please ensure your operating system has all prerequisite packages installed and is compatible with R GMMAT and Matrix packages. Second, please convert your genotype files into .gds format before running the local ancestry analysis. The function processes variants in batches (default = 100 SNPs per batch) to account for local ancestry effects. It automatically handles missing values in the input data and returns a dataframe containing minor allele frequencies (MAF), p-values, and other relevant statistics.

```
### Establish a null model via glmmkin()
pheno_file_path<-system.file('estdata','pheno',package='LOCATOR')
null.model<-GMMAT::glmmkin(y~age+sex,data=pheno,id='IID')
### Run LOCATOR
gds<-system.file('estdata','geno.gds',package='LOCATOR')
lac_file_path<-system.file('estdata','LACs_oro',package='LOCATOR')
get.la.aware.score(null.obj=null.model,geno.file=gds,outfile=%placeholder%,LAC.file=lac_file_path,id.type='character')
```

### Step 4: Map phenotypic distributions onto global PCs
The plot.features() function visualizes categorical phenotypic distributions on the global PC canvas, helping illustrate relationships between phenotypes, global PCs, and LACs. It generates pie charts representing phenotypic distributions in terms of local ancestry patterns at specific positions.
We build plot.features() on the scatterpie and ggplot2 R packages, so please require these packages in advance. It uses shared coordinates with global PCs and LACs. Since these coordinates are typically more concentrated than phenotypic distributions, please use the y.jitter argument to adjust pie chart radius and if.jitter to modify chart positions to prevent overlapping. Example code for generating phenotypic distribution plots is provided for further illustration and customization.

```
### The genetic position at 12:48339919 with the minimum p-value (p = 0.0002) would be plotted
### Extract corresponding LACs information
lac_file_path<-system.file('estdata','LACs',package='LOCATOR')
con<-file(lac_file_path,’rb’)
readBin(con=con,what=integer(),n=1,size=4)  # number of PCs: 5
readBin(con=con,what=integer(),n=1,size=4)  # number of samples: 200
readBin(con=con,what=character(),n=200,size=4)  # sample IDs
readBin(con=con, what=integer(),n=1,size=4)  # number of breakpoints: 100
pos<-readBin(con=con,what=integer(),n=100,size=4)
which.max(pos<=48336619)  ## identify the corresponding tract
la_coord<-readBin(con,what=numeric(),n=200*5,size=4)
la_coord<-as.data.frame.matrix(matrix(la_coord,ncol=5))
colnames(la_coord)<-paste0('PC',1:5,'_la')
close(con=con)
### Build the meta dataframe
plot_file<-cbind(gc_file[,1:3],la_coord[,1:2])
plot_file<-mutate(plot_file,anc=ga_file$group,pheno=pheno$y)
plot.features(input.file=plot_file,PC1=PC1,PC2=PC2,anc=anc,PC1_la=PC1_la,PC2_la=PC2_la,pheno=y,pheno.name='pheno',phenotype.name='phenotype',y.scale=4000,r.x.axis=0,r.y.axis=0.2,r.n.level=3)
```

### Caveats and cautions:
1. Please be aware that current local ancestry inference software provides universal breakpoints across all tested samples, rather than individual-specific breakpoints. While these tools generally achieve high accuracy, some errors may occur in complex genomic regions, particularly those with high linkage disequilibrium (LD). Please be careful when dealing with SNPs in these specific regions.
2. We currently set byte size = 4 for integers and numerics storage in binary files, with a maximum value of 2^32 - 1 = 4,294,967,295. For chromosomes with genetic positions exceeding 4.3 billion, please subset your genotype files accordingly. Meanwhile, please avoid combining genetic positions from different chromosomes in a single file to prevent position index interference. If you use double numbers for sample labeling, please switch them into either characters, integers, or numerics that meet with the 4-byte-size storage requirement.
Please check out your local ancestry input files in advance. Some local ancestry inference outputs may contain certain missing values, especially with multiple reference populations (e.g. more than 10 reference ancestries). The get.la.aware.score() will fill missing values with prior local ancestry values. However, the process will halt if no prior local ancestry/LACs are available or if missing values exceed 10%.
3. The toolkit currently supports haplotype-based association tests for autosomes only. Please ensure your genotype input files are properly converted to haplotypic format. Please contact us with specific file format requirements to extend haplotype association tests beyond autosomes. Thank you.
4. Please be aware that binary file compatibility is platform-dependent. The endian (byte-order) argument is not specified by default when reading/writing binary files to accommodate various operating platforms. To maintain consistency, please perform local ancestry analysis on the same platform, or specify endian when necessary. For more details, visit https://www.rdocumentation.org/packages/base/versions/3.6.2/topics/readBin.
