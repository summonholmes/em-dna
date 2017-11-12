def input_start_iter():
	input_start_iter_ = int(input("Use 500 iterations to perform the E-M steps? ('1' for yes or '0' for no): "))
	number_iter = 0
	if input_start_iter_ == 1:
		number_iter = 500
	elif input_start_iter_ == 0:
		user_input_rand = int(input("Please specify a value for the number of iterations: "))
		number_iter = user_input_rand
	else:
		print("Unknown option.  Now exiting")
		exit(0)
	return number_iter
