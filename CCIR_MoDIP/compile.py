import csv
import numpy as np
for mth in range(1,13):
    with open('ccir' + str(10 + mth) + '.txt') as f:
        data = []
        for row in csv.reader(f, delimiter=' '):
            row = [num for num in row if num != '']  # filter
            data = data + row
        assert (len(data) == 2858)

        F2data = data[:1976]
        F2 = np.reshape(np.array(F2data), (2, 76, 13))

        Fm3data = data[1976:]
        Fm3 = np.reshape(np.array(Fm3data), (2, 49, 9))

        F2list = F2.tolist()
        Fm3list = Fm3.tolist()

        with open('F2.py', 'a') as F2out:
            F2out.writelines('CCIR' + str(10 + mth) + '_F2 = [\n' )
            for matrix in F2list:
                F2out.writelines('[\n')
                F2out.writelines('[' + ','.join(row) + '],\n' for row in matrix)
                F2out.writelines('],\n')
            F2out.writelines(']\n\n' )

        with open('Fm3.py', 'a') as Fm3out:
            Fm3out.writelines('CCIR' + str(10 + mth) + '_Fm3 = [\n' )
            for matrix in Fm3list:
                Fm3out.writelines('[\n')
                Fm3out.writelines('[' + ','.join(row) + '],\n' for row in matrix)
                Fm3out.writelines('],\n')
            Fm3out.writelines(']\n\n' )

with open('F2.py', 'a') as F2out:
    F2out.writelines('F2 = {1: CCIR11_F2, 2: CCIR12_F2, 3: CCIR13_F2, 4: CCIR14_F2, 5: CCIR15_F2, 6: CCIR16_F2, 7: CCIR17_F2, 8: CCIR18_F2, 9: CCIR19_F2, 10: CCIR20_F2, 11: CCIR21_F2, 12: CCIR22_F2}')

with open('Fm3.py', 'a') as Fm3out:
    Fm3out.writelines('Fm3 = {1: CCIR11_FM3, 2: CCIR12_FM3, 3: CCIR13_FM3, 4: CCIR14_FM3, 5: CCIR15_FM3, 6: CCIR16_FM3, 7: CCIR17_FM3, 8: CCIR18_FM3, 9: CCIR19_FM3, 10: CCIR20_FM3, 11: CCIR21_FM3, 12: CCIR22_FM3}')

