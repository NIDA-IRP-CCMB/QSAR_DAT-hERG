import getopt
import os
import sys
import time


def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hm:x:", ["model=", "method="])
    except getopt.GetoptError:
        print(sys.argv[0] + ' -m <model> -x <method>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0] + ' -c <class>  -m <model>')
            sys.exit()
        elif opt in ('-m', '--model'):
            mode = arg
        elif opt in ('-x', '--method'):
            ml = arg

    return mode, ml


if __name__ == "__main__":
    mode, ml = main(sys.argv[1:])
    print("running build model %s %s" % (mode,ml))

    if mode.startswith('reg'):
        mode = 'reg'
        os.system("time python ${COREPATH}/run_buildmodel.py -m %s -x %s -n 1 -r 0 " % (mode, ml) +
                  "-i dataset -b pubdata_40 -t 0.00 -o output_%s_%s \n"
                  % (mode, ml))
    elif mode.startswith('class'):
        mode = 'class'
        os.system("time python ${COREPATH}/run_buildmodel.py -m %s -x %s -n 1 -r 0" % (mode, ml) +
                  "-i dataset -b pubdata_class_40 -t 0.00 -o output_%s_%s  \n"
                  % (mode, ml))
    print("running prediction %s %s" % (mode, ml))
    os.system("\n time python ${COREPATH}/run_prediction.py -m %s -x %s -r 0 -n 1 " % (mode, ml) +
              "-i output_%s_%s -b testset_pairs -d validation_dataset -t 0.00 -o output"
              % (mode, ml))
