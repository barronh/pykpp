

if __name__ == '__main__':
    from optparse import OptionParser
    from test import testit
    parser = OptionParser()
    parser.set_usage("""Usage: python -m pykpp.test [-v] model1,model2""")

    parser.add_option("-v", "--verbose", dest = "verbose", action = "store_true", default = False, help = "Show extended output")
    options, args = parser.parse_args()
    testit(*args, verbose = options.verbose)