#!/usr/bin/env bash
#
# Pipeline: BUSCO completeness assessment for candidate tardigrade
# reference genome assemblies.
#
# Prerequisites (install once):
#   1) Miniconda/Mamba installed
#   2) NCBI Datasets CLI (to download the genomes)
#   3) BUSCO v5.x (via conda-forge/bioconda)
#
# Suggested installation (run once, BEFORE running this script):
#
#   conda create -n busco_env -c conda-forge -c bioconda busco=5.7.1 -y
#   conda create -n ncbi_datasets -c conda-forge ncbi-datasets-cli -y
#
# This script activates each conda environment automatically when it
# needs it (you don't need to activate anything manually). You just
# need conda installed and both environments already created.
#
# ---------------------------------------------------------------------

set -euo pipefail

# Conda environment names (adjust if you used different names)
NCBI_ENV="ncbi_datasets"
BUSCO_ENV="busco_env"

# Load conda activate/deactivate functions inside the script.
# Without this, "conda activate" fails in non-interactive scripts.
CONDA_BASE="$(conda info --base)"
source "${CONDA_BASE}/etc/profile.d/conda.sh"

### 0. CONFIGURATION ####################################################

# Base working directory
BASE_DIR="$(pwd)/tardigrade_busco"
GENOMES_DIR="${BASE_DIR}/genomes"
RESULTS_DIR="${BASE_DIR}/busco_results"

# BUSCO lineage to use (same for all 5 species -> comparable results)
# CONFIRMED: the NCBI report for Pam. metropolitanus used
# metazoa_odb10 (n=954), with BUSCO v4.1.4. This pipeline runs BUSCO
# v5.7.1 for all 5 species -- INCLUDING Pam. metropolitanus again --
# so that all 5 values come from the same version and are directly
# comparable to each other. The original NCBI value (C:85.0%[S:81.3%,
# D:3.7%],F:1.8%,M:13.2%,n:954) is kept only as a historical reference,
# it is not used in the final comparative table.
LINEAGE="metazoa_odb10"

# Number of cores to use
THREADS=8

# BUSCO mode
MODE="genome"

# NCBI accessions (GCF/GCA) for the 5 candidate references
declare -A GENOMES=(
  ["Pam_metropolitanus"]="GCF_019649055.1"
  ["Hyp_dujardini"]="GCA_001579985.1"
  ["Ram_varieornatus"]="GCA_001949185.1"
  ["Hyp_exemplaris"]="GCA_002082055.1"
  ["Par_richtersi"]="GCA_034698045.1"
)

mkdir -p "${GENOMES_DIR}" "${RESULTS_DIR}"

### 1. DOWNLOAD GENOMES (NCBI Datasets) #################################

echo ">>> Activating '${NCBI_ENV}' environment for download..."
conda activate "${NCBI_ENV}"

echo ">>> Downloading genomes from NCBI Datasets..."

for species in "${!GENOMES[@]}"; do
  accession="${GENOMES[$species]}"
  species_dir="${GENOMES_DIR}/${species}"

  if [[ -f "${species_dir}/${species}.fasta" ]]; then
    echo "  - ${species} (${accession}) already downloaded, skipping."
    continue
  fi

  echo "  - Downloading ${species} (${accession})..."
  mkdir -p "${species_dir}"

  datasets download genome accession "${accession}" \
    --include genome \
    --filename "${species_dir}/${accession}.zip"

  unzip -o -q "${species_dir}/${accession}.zip" -d "${species_dir}/unzipped"

  # The FASTA ends up at ncbi_dataset/data/<accession>/*.fna
  fasta_path=$(find "${species_dir}/unzipped" -name "*.fna" | head -n 1)
  cp "${fasta_path}" "${species_dir}/${species}.fasta"

  # Clean up intermediate files
  rm -rf "${species_dir}/unzipped" "${species_dir}/${accession}.zip"
done

echo ">>> Download complete."
conda deactivate
echo

### 2. RUN BUSCO FOR EACH GENOME ########################################

echo ">>> Activating '${BUSCO_ENV}' environment for BUSCO..."
conda activate "${BUSCO_ENV}"

# FIXED gene predictor for all 5 species. Metaeuk was chosen instead of
# the default (Miniprot) because Miniprot showed low prediction identity
# and internal stop codons for Pam_metropolitanus. Using the same
# predictor across all 5 species is essential for BUSCO scores to be
# directly comparable to each other -- mixing predictors depending on
# which one failed would introduce a bigger inconsistency than the one
# being fixed.
PREDICTOR_FLAG="--metaeuk"
PREDICTOR_NAME="metaeuk"

echo ">>> Running BUSCO (lineage: ${LINEAGE}, predictor: ${PREDICTOR_NAME})..."

# IMPORTANT: the previous Pam_metropolitanus result (run with Miniprot)
# is moved aside so all 5 species end up under the same predictor. The
# already-downloaded genome is NOT re-downloaded, only the BUSCO run is
# repeated for that species.
if [[ -d "${RESULTS_DIR}/busco_Pam_metropolitanus" ]]; then
  echo "  - Found a previous Pam_metropolitanus result (Miniprot)."
  echo "    Moving it to *_miniprot_backup to keep it, then re-running with Metaeuk."
  mv "${RESULTS_DIR}/busco_Pam_metropolitanus" \
     "${RESULTS_DIR}/busco_Pam_metropolitanus_miniprot_backup"
fi

declare -A PREDICTOR_USED=()

for species in "${!GENOMES[@]}"; do
  fasta="${GENOMES_DIR}/${species}/${species}.fasta"
  out_name="busco_${species}"

  if [[ -d "${RESULTS_DIR}/${out_name}" ]]; then
    echo "  - ${species} already has BUSCO results (${PREDICTOR_NAME}), skipping."
    PREDICTOR_USED["${species}"]="${PREDICTOR_NAME} (already_existing)"
    continue
  fi

  echo "  - Running BUSCO on ${species} (${PREDICTOR_NAME})..."

  # Temporarily disable "exit on error": if BUSCO fails on one species,
  # we don't want it to abort the rest -- it gets logged as failed and
  # the loop moves on to the next species instead.
  set +e
  busco \
    -i "${fasta}" \
    -o "${out_name}" \
    -l "${LINEAGE}" \
    -m "${MODE}" \
    -c "${THREADS}" \
    --out_path "${RESULTS_DIR}" \
    ${PREDICTOR_FLAG} \
    -f
  status=$?
  set -e

  if [[ ${status} -eq 0 ]]; then
    PREDICTOR_USED["${species}"]="${PREDICTOR_NAME}"
  else
    echo "  !! ${species} failed with ${PREDICTOR_NAME} (exit ${status})."
    echo "     Check manually: ${RESULTS_DIR}/${out_name}/logs/busco.log"
    PREDICTOR_USED["${species}"]="FAILED"
  fi
done

echo ">>> BUSCO finished for all species (see details above)."
echo ">>> Predictor used per species:"
for species in "${!PREDICTOR_USED[@]}"; do
  echo "    - ${species}: ${PREDICTOR_USED[$species]}"
done
conda deactivate
echo

### 3. CONSOLIDATED SUMMARY #############################################
# (needs no special environment, just standard grep/awk)

SUMMARY_FILE="${RESULTS_DIR}/busco_summary_all_species.tsv"
echo -e "species\taccession\tcomplete_single\tcomplete_duplicated\tfragmented\tmissing\tn_genes_total\tpredictor" > "${SUMMARY_FILE}"

for species in "${!GENOMES[@]}"; do
  accession="${GENOMES[$species]}"
  predictor="${PREDICTOR_USED[$species]:-unknown}"
  short_summary=$(find "${RESULTS_DIR}/busco_${species}" -name "short_summary*.txt" | head -n 1)

  if [[ -z "${short_summary}" ]]; then
    echo "  ! No summary found for ${species}, check log."
    echo -e "${species}\t${accession}\tNA\tNA\tNA\tNA\tNA\t${predictor}" >> "${SUMMARY_FILE}"
    continue
  fi

  # Extract values from BUSCO's short_summary.txt
  single=$(grep "Complete and single-copy" "${short_summary}" | awk '{print $1}')
  dup=$(grep "Complete and duplicated"   "${short_summary}" | awk '{print $1}')
  frag=$(grep "Fragmented"               "${short_summary}" | awk '{print $1}')
  miss=$(grep "Missing"                  "${short_summary}" | awk '{print $1}')
  total=$(grep "Total BUSCO groups"      "${short_summary}" | awk '{print $1}')

  echo -e "${species}\t${accession}\t${single}\t${dup}\t${frag}\t${miss}\t${total}\t${predictor}" >> "${SUMMARY_FILE}"
done

echo ">>> Consolidated summary written to: ${SUMMARY_FILE}"
echo ">>> Check that file to complete your LaTeX table."

# Final warning if any species failed to complete with either predictor
failed_any=0
for species in "${!PREDICTOR_USED[@]}"; do
  if [[ "${PREDICTOR_USED[$species]}" == "FAILED" ]]; then
    failed_any=1
    echo "!!! WARNING: ${species} could not complete even with ${PREDICTOR_NAME}."
  fi
done

if [[ ${failed_any} -eq 1 ]]; then
  echo "!!! Please check the species marked as FAILED above manually."
fi