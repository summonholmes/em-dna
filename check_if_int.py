def check_if_int(input_start_align):
    try:
        input_start_align = int(input_start_align)
        return input_start_align
    except:
        print("You didn't provide an integer ya turkey!")
        exit(0)