from NequickG import GalileoBroadcast, NEQTime
from NequickG_global import NequickG_global
import csv
import matplotlib.pyplot as plt
import numpy as np

def run(table_type):
    """
    For each entry in validation table, compute sTEC
    Validation run takes a few minutes to complete so results are saved to file.
    :param table_type: (str)
    :return:
    """
    print "Validating against", table_type, "Solar Activity data table"

    # Hard coded solar indices
    if table_type == "High":
        BX = GalileoBroadcast(236.831641, -0.39362878, 0.00402826613)
    elif table_type == "Medium":
        BX = GalileoBroadcast(121.129893, 0.351254133, 0.0134635348)
    elif table_type == "Low":
        BX = GalileoBroadcast(2.580271, 0.127628236, 0.0252748384)
    else:
        raise ValueError('table_type argument must be either "High", "Medium" or "Low"')

    # output variables
    sTECs_expected = []
    sTECs_computed = []

    with open('./Validation/' + table_type + '_reference.dat') as infile:
        with file('./Validation/' + table_type + '_output.dat', 'w') as outfile:
            writer = csv.writer(outfile, delimiter = ' ')
            reader = csv.reader(infile, delimiter = ' ')
            for row in reader:
                row = map(float, row)

                mth = int(row[0])
                UT = int(row[1])

                lon1 = row[2]
                lat1 = row[3]
                h1 = row[4] / 1000 # change m to km

                lon2 = row[5]
                lat2 = row[6]
                h2 = row[7] / 1000 # change m to km

                time = NEQTime(mth, UT)
                NEQ_global = NequickG_global(time, BX)
                stec = NEQ_global.sTEC(h1, lat1, lon1, h2, lat2, lon2) / 10**16

                sTECs_computed.append(stec)
                sTECs_expected.append(row[8])
                row.append (stec)

                writer.writerow(row)

                print "Input: ", row[:7]
                print "---Expected: ", row[8], "---Obtained: ", row[9]
    return sTECs_expected, sTECs_computed

def get_computed(table_type):
    try:
        df = np.loadtxt('Validation/' + table_type + '_output.dat', delimiter=' ')
    except IOError:
        raise ValueError('Does .out file exist? table_type argument must be either "High", "Medium" or "Low"')

    sTECs_expected, sTECs_computed = df[:,8], df[:,9]

    return sTECs_expected, sTECs_computed


def compare(sTECs_expected, sTECs_computed):
    # Change List type to np.ndarray type
    sTECs_expected = np.array(sTECs_expected)
    sTECs_computed = np.array(sTECs_computed)

    # analysis
    residual = sTECs_expected - sTECs_computed
    A = np.vstack([sTECs_expected, np.ones(len(sTECs_expected))]).T
    (m,c) = np.linalg.lstsq(A,sTECs_computed)[0]

    plt.subplot(2,1,1)
    plt.plot(sTECs_expected, sTECs_computed, 'ro')
    plt.plot(sTECs_expected, m * sTECs_expected + c)
    # plt.xlabel('Expected sTEC')
    plt.ylabel('Computed sTEC')
    plt.grid()

    plt.subplot(2,1,2)
    plt.plot(sTECs_expected, residual, 'ro')
    plt.xlabel('Expected sTEC')
    plt.ylabel('Residual')
    plt.grid()
    plt.show()

# sTECs_expected, sTECs_computed = run('Medium')
sTECs_expected, sTECs_computed = get_computed('Medium')
compare(sTECs_expected, sTECs_computed)
