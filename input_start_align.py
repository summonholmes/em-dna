from check_if_int import check_if_int


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