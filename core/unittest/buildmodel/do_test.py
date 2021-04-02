import getopt
import os
import sys
import time

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hc:m:", ["class=", "model="])
    except getopt.GetoptError:
        print(sys.argv[0] + ' -c <class> -m <model>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0] + ' -c <class>  -m <model>')
            sys.exit()
        elif opt in ('-c', '--class'):
            c = arg
        elif opt in ('-m', '--model'):
            ml = arg

    return c, ml


if __name__ == "__main__":
    c, ml = main(sys.argv[1:])

    if c.startswith('reg'):
        c = 'reg'
    elif c.startswith('class'):
        c = 'class'
    print("running build model " + c + " " + ml)
    os.system("time python core/run_buildmodel.py -c %s -m %s " % (c, ml) +
              "-i dataset_rDAT_inhibitor_Ki -b pubdata -o output_%s_%s -t 0.00 -n 1 \n"
              % (c, ml))
    print("running prediction " + c + " " + ml)
    os.system("\n time python core/run_prediction.py -c %s " % c +
              "-i output_%s_%s -o output -a -d validation_dataset -b testset_pairs"
              % (c, ml))
