#!/bin/zsh
#v1.1.0
show_usage() {
    echo "Usage: ./run_ngs.sh [-n <Path to the reference sequence>] <File with samples>"
    echo "Example: ./run_ngs.sh -n ref_seq.fasta samples.txt"
}

error_exit() {
    echo "Error: $1" >&2
    exit "${2:-1}"
}

while getopts "n:" opt; do
    case $opt in
        n) path_to_nat_sequence="$OPTARG";;
        \?) echo "Invalid argument: -$OPTARG" >&2; show_usage; exit 1;;
        :) echo "argument -$OPTARG expects defintion." >&2; show_usage; exit 1;;
    esac
done

# Move the position of the argument to take the "optind" into account -> $1 = samples
shift $((OPTIND - 1))
sample_file="$1"

# Checks whether -n ist defined
if [ -z "$path_to_nat_sequence" ]; then
    show_usage
    error_exit "The options -n and -f must be defined." 1
fi

# Checks whether the number of arguments is correct
if [ "$#" -ne 1 ]; then
    show_usage
    error_exit "Incorrect number of arguments. Please enter frameshift, path to the native sequence and sample file." 1
fi

# Checks whether the path_to_nat_sequence file exists
if [ -z "$path_to_nat_sequence" ] || [ ! -f "$path_to_nat_sequence" ]; then
    show_usage
    error_exit "The reference sequence file '$path_to_nat_sequence' does not exist." 1
fi

# Checks whether all sample_names are unique
while IFS=',' read -r line || [[ -n "$line" ]]; do
    IFS=',' read -r -A fields <<< "$line"
    sample_name="${fields[-5]}"
    if [[ " ${processed_samples[@]} " =~ " $sample_name " ]]; then
        error_exit "Duplicate sample_name '$sample_name' found." 1
    else
        processed_samples+=("$sample_name")
    fi
done < "$sample_file"

mkdir -p ./demultiplex || error_exit "Demultiplex directory could not be created."
mkdir -p ./results || error_exit "Results directory could not be created."

demultiplex_sequence() {
    # Function that does the demultiplexing
    local fasta_file=$1
    local sample_name=$2
    local bc_start=$3
    local bc_end=$4
    local length=$5
    local output_seq="demultiplex/${sample_name}.seq"
    local output_fasta="demultiplex/${sample_name}.fasta"

    #demultiplex
    grep -o "$bc_start.*$bc_end" "$fasta_file.fasta" | sed -e "s/.*$bc_start\(.*\)$bc_end.*/\1/" | awk "length(\$0) ==$length" > "$output_seq" || error_exit "Demultiplex error for ${sample_name}."
    
    #.seq to .fasta
    if [[ -e "$output_fasta" ]]; then
        error_exit "The file '$output_fasta' already exists." 2
    else
        counter=0
        while IFS= read -r line; do
            let counter++
            echo ">${sample_name}_${counter}" >> "$output_fasta"
            echo "$line" >> "$output_fasta"
        done < "$output_seq"
        seq_count=$(wc -l < "$output_fasta")
        echo "${sample_name} demultiplexed - Number of sequences: ${seq_count}"
    fi
}

calculate_reading_frame() {
    # Function takes the first sequence line, divide by 3 using a modulo operation and returns the the ReadingFrame.
    local sample_name=$1
    
    length=$(sed -n '2p' "demultiplex/${sample_name}_S1.fasta" | tr -d '\n' | wc -c)

    remainder=$((length % 3))
    case $remainder in
        0)
            reading_frame=1
            ;;
        1)
            reading_frame=2
            ;;
        2)
            reading_frame=3
            ;;
    esac
}

echo "Start Demultiplexing..."

while IFS=',' read -r line || [[ -n "$line" ]]; do
    plotting_successful=true
    ((line_number++))
    adjusted_line_number=$((line_number - 1))

    IFS=',' read -rA fields <<< "$line"
    num_vars=${#fields[@]}
    num_fasta=$((num_vars - 5))

    for ((i = 1; i <= num_fasta; i++)); do
        demultiplex_sequence $fields[$i] $fields[-5]_S$i $fields[-4] $fields[-3] $fields[-2]
    done

    echo "Plotting" $fields[-5]
    calculate_reading_frame $fields[-5] #returns $reading_frame
    echo $reading_frame
    python plot.py "$reading_frame" "$path_to_nat_sequence" "$sample_file" "$adjusted_line_number" || { error_exit "Python script error."; ; plotting_successful=false; }
done < "$sample_file"

if [ "$plotting_successful" = true ]; then
    echo "Plotting done!"
fi

# __author__ = "Tom U. Schlegel"
# __contact__ = "tom.schlegel@uni-leipzig.de"
# __license__ = "GNU GPLv3"
# __version__ = "1.1"