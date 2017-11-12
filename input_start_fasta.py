from Bio import SeqIO


def input_start_fasta():
	user_choice = 1
	# fasta_file_id = []
	fasta_file_seq = []
	while user_choice != 0:
		fasta_read = int(input("\nImport fasta file? ('1' for yes or '0' for no): "))
		if fasta_read == 1:
			while True:
				user_fasta = input("Please specify the path of the fasta file, or press '0' to continue: ")
				if user_fasta == '0':
					break
				else:
					try:
						for record in SeqIO.parse(user_fasta, "fasta"):
							# fasta_file_id.append(str(record.id))
							fasta_file_seq.append(str(record.seq))
					except IOError:
						print("Unable to open file.  Try again.")
		elif fasta_read == 0:
			break
		else:
			print("Unknown option, now exiting.")
			exit(0)
	print("\nNow preparing for E-M...")
	print("M = Max Motif\nS = Max Score\nP = Max Position\n# = Max Sequence #\nSS = Max Sum Scores")
	print("\n\t\t\t\t\t\tM:\t\t S:\t\t  P:  #:\tSS:")
	return fasta_file_seq
