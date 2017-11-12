def input_start_align():
	input_start_align_ = int(input("Use 50 initial random starting alignments? ('1' for yes or '0' for no): "))
	random_start_align = 0
	if input_start_align_ == 1:
		random_start_align = 50
	elif input_start_align_ == 0:
		user_input_rand = int(input("Please specify a value for the initial random starting alignments: "))
		random_start_align = user_input_rand
	else:
		print("Unknown option.  Now exiting")
		exit(0)
	return random_start_align
