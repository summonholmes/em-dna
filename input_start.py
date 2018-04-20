from check_if_exception import check_if_int, check_if_fasta


def input_start_align():
    input_start_align = input(
        "Use 50 initial random starting alignments? ('1' for yes or '0' for no): "
    )
    input_start_align = check_if_int(input_start_align)
    if input_start_align == 1:
        random_start_align = 50
    elif input_start_align == 0:
        user_input_rand = input(
            "Please specify a value for the initial random starting alignments: "
        )
        user_input_rand = check_if_int(user_input_rand)
        random_start_align = user_input_rand
    else:
        print("Unknown option.  Now exiting")
        exit(0)
    return random_start_align


def input_start_iter():
    input_start_iter = input(
        "Use 500 iterations to perform the E-M steps? ('1' for yes or '0' for no): "
    )
    input_start_iter = check_if_int(input_start_iter)
    if input_start_iter == 1:
        number_iter = 500
    elif input_start_iter == 0:
        number_iter = input(
            "Please specify a value for the number of iterations: ")
        number_iter = check_if_int(number_iter)
    else:
        print("Unknown option.  Now exiting")
        exit(0)
    return number_iter


def input_start_fasta():
    fasta_file_seq = []
    user_fasta_path = input("Please specify the path of the fasta file: ")
    fasta_file_seq = check_if_fasta(fasta_file_seq, user_fasta_path)
    return fasta_file_seq