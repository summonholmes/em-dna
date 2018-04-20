from check_if_int import check_if_int


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