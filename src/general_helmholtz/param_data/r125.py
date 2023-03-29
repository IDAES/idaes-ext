from helmholtz_parameters import WriteParameters


def main():
    we = WriteParameters(parameters="r125.json")
    we.write()


if __name__ == "__main__":
    main()
