library(Peptides)
library(tools)

get_seq_fragment <- function(seq, idx) {
	start_pos <- max(1, idx - 3)
	end_pos <- min(nchar(seq), idx + 3)

	if ((end_pos - start_pos + 1) < 7) {
		if (start_pos == 1) {
			end_pos <- min(nchar(seq), start_pos + 6)
		}
		else if (end_pos == nchar(seq)) {
			start_pos <- max(1, end_pos - 6)
		}
	}

	return(substring(seq, start_pos, end_pos))
}


results_df <- data.frame(
	pdbid = character(),
	chain = character(),
	labeling_pos = integer(),
	KF1 = numeric(),
	KF2 = numeric(),
	KF3 = numeric(),
	KF4 = numeric(),
	KF5 = numeric(),
	KF6 = numeric(),
	KF7 = numeric(),
	KF8 = numeric(),
	KF9 = numeric(),
	KF10 = numeric(),
	PP1 = numeric(),
	PP2 = numeric(),
	PP3 = numeric(),
	VSHE1 = numeric(),
    VSHE2 = numeric(),
    VSHE3 = numeric(),
    VSHE4 = numeric(),
    VSHE5 = numeric(),
    VSHE6 = numeric(),
    VSHE7 = numeric(),
    VSHE8 = numeric(),
	F1 = numeric(),
    F2 = numeric(),
    F3 = numeric(),
    F4 = numeric(),
    F5 = numeric(),
    F6 = numeric(),
	Z1 = numeric(),
    Z2 = numeric(),
    Z3 = numeric(),
    Z4 = numeric(),
    Z5 = numeric(),
    T1 = numeric(),
    T2 = numeric(),
    T3 = numeric(),
    T4 = numeric(),
    T5 = numeric(),
    ST1 = numeric(),
    ST2 = numeric(),
    ST3 = numeric(),
    ST4 = numeric(),
    ST5 = numeric(),
    ST6 = numeric(),
    ST7 = numeric(),
    ST8 = numeric(),
    FP1 = numeric(),
    FP2 = numeric(),
    FP3 = numeric(),
    FP4 = numeric(),
    FP5 = numeric(),
    FP6 = numeric(),
    FP7 = numeric(),
    FP8 = numeric(),
	BL1 = numeric(),
    BL2 = numeric(),
    BL3 = numeric(),
    BL4 = numeric(),
    BL5 = numeric(),
    BL6 = numeric(),
    BL7 = numeric(),
    BL8 = numeric(),
	BL9 = numeric(),
    BL10 = numeric(),

	stringsAsFactors = FALSE
)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Please provide PDB ID and chain as arguments.")
}

pdbid <- args[1]
chain <- args[2]
fasta <- paste0('./', pdbid, '_', chain, '.fasta')
fasta_contents <- readLines(fasta)
fasta_seq <- fasta_contents[2]

hrpf_file <- paste0('./', pdbid, '_', chain, '_HRPF.csv')
lines <- readLines(hrpf_file)
num_lines <- length(lines)
label_position_vector <- integer(num_lines)

for (i in seq_along(lines)){
	split_line <- strsplit(lines[i], "\t")[[1]]
	labeling_pos <- as.integer(split_line[1])
	fasta_idx <- labeling_pos
	print(fasta_idx)
	label_position_vector[i] <- as.integer(fasta_idx)
}
for (i in seq_along(label_position_vector)) {
	seq_frag <- get_seq_fragment(fasta_seq, label_position_vector[i])
	kF_results <- kideraFactors(seq_frag)
	kF_results_vector <- unlist(kF_results)
	
	cP_results <- crucianiProperties(seq_frag)
	cP_results_vector <- unlist(cP_results)

	vhse_results <- vhseScales(seq_frag)
	vhse_results_vector <- unlist(vhse_results)

	fV_results <- fasgaiVectors(seq_frag)
	fV_results_vector <- unlist(fV_results)

	z_results <- zScales(seq_frag)
	z_results_vector <- unlist(z_results)

	t_results <- tScales(seq_frag)
	t_results_vector <- unlist(t_results)

	st_results <- stScales(seq_frag)
	st_results_vector <- unlist(st_results)

	fP_results <- protFP(seq_frag)
	fP_results_vector <- unlist(fP_results)

	blosum_results <- blosumIndices(seq_frag)
	blosum_results_vector <- unlist(blosum_results)

	row <- data.frame(
		pdbid = pdbid,
		chain = chain,
		labeling_pos = label_position_vector[i],
		KF1 = kF_results_vector["KF1"],
		KF2 = kF_results_vector["KF2"],
		KF3 = kF_results_vector["KF3"],
		KF4 = kF_results_vector["KF4"],
		KF5 = kF_results_vector["KF5"],
		KF6 = kF_results_vector["KF6"],
		KF7 = kF_results_vector["KF7"],
		KF8 = kF_results_vector["KF8"],
		KF9 = kF_results_vector["KF9"],
		KF10 = kF_results_vector["KF10"],
		PP1 = cP_results_vector["PP1"],
		PP2 = cP_results_vector["PP2"],
		PP3 = cP_results_vector["PP3"],
		VHSE1 = vhse_results_vector['VHSE1'],
		VHSE2 = vhse_results_vector['VHSE2'],
		VHSE3 = vhse_results_vector['VHSE3'],
		VHSE4 = vhse_results_vector['VHSE4'],
		VHSE5 = vhse_results_vector['VHSE5'],
		VHSE6 = vhse_results_vector['VHSE6'],
		VHSE7 = vhse_results_vector['VHSE7'],
		VHSE8 = vhse_results_vector['VHSE8'],
		F1 = fV_results_vector['F1'],
		F2 = fV_results_vector['F2'],
		F3 = fV_results_vector['F3'],
		F4 = fV_results_vector['F4'],
		F5 = fV_results_vector['F5'],
		F6 = fV_results_vector['F6'],
		Z1 = z_results_vector['Z1'],
		Z2 = z_results_vector['Z2'],
		Z3 = z_results_vector['Z3'],
		Z4 = z_results_vector['Z4'],
		Z5 = z_results_vector['Z5'],
		T1 = t_results_vector['T1'],
		T2 = t_results_vector['T2'],
		T3 = t_results_vector['T3'],
		T4 = t_results_vector['T4'],
		T5 = t_results_vector['T5'],
		ST1 = st_results_vector['ST1'],
		ST2 = st_results_vector['ST2'],
		ST3 = st_results_vector['ST3'],
		ST4 = st_results_vector['ST4'],
		ST5 = st_results_vector['ST5'],
		ST6 = st_results_vector['ST6'],
		ST7 = st_results_vector['ST7'],
		ST8 = st_results_vector['ST8'],
		FP1 = fP_results_vector['ProtFP1'],
		FP2 = fP_results_vector['ProtFP2'],
		FP3 = fP_results_vector['ProtFP3'],
		FP4 = fP_results_vector['ProtFP4'],
		FP5 = fP_results_vector['ProtFP5'],
		FP6 = fP_results_vector['ProtFP6'],
		FP7 = fP_results_vector['ProtFP7'],
		FP8 = fP_results_vector['ProtFP8'],
		BL1 = blosum_results_vector['BLOSUM1'],
		BL2 = blosum_results_vector['BLOSUM2'],
		BL3 = blosum_results_vector['BLOSUM3'],
		BL4 = blosum_results_vector['BLOSUM4'],
		BL5 = blosum_results_vector['BLOSUM5'],
		BL6 = blosum_results_vector['BLOSUM6'],
		BL7 = blosum_results_vector['BLOSUM7'],
		BL8 = blosum_results_vector['BLOSUM8'],
		BL9 = blosum_results_vector['BLOSUM9'],
		BL10 = blosum_results_vector['BLOSUM10'],
		stringsAsFactors = FALSE
	)
		results_df <- rbind(results_df, row)
}
write.csv(results_df, paste0(pdbid, "_", chain, "_R_features.csv"), row.names = FALSE)
