
python3 wgetter.py

var="$(cat temp.txt)"

if [ ! -z "$var" ];
then
	wget "$var"
	gunzip *.gz
	mv *.faa OrnAna.fa
	makeblastdb -in OrnAna.fa -dbtype prot -parse_seqids
	blastx -query T1.fa -db OrnAna.fa -outfmt 6 -out total.txt
	sort -k1,1 -k11,11g total.txt | sort --merge -u  -k1,1 > total_sort.txt
	cat total_sort.txt | awk '{print $2}' > ensemble_ids.txt

	python3 script_ncbi.py ensemble_ids.txt

	rm O*
	rm tot*
	rm temp.txt
	rm esummary*
	rm ensemble_ids.txt
fi
