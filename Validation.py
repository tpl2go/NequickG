from NequickG import NequickG, NequickG_global, Position, GalileoBroadcast, NEQTime
import csv
# High Solar Activity

def highsolaractivity():
    BX = GalileoBroadcast(236.831641, -0.39362878, 0.00402826613)
    # Month, Time, Station longitude, Station latitude, Station height, Statellite longitude, statlite latitude, satellite height
    with open('./Validation/High_reference.dat') as infile:
        with open('./Validation/High_output.dat', 'w') as outfile:
            writer = csv.writer(outfile, delimiter = ' ')
            reader = csv.reader(infile, delimiter = ' ')
            for row in reader:
                row = map(float, row)
                print row

                mth = int(row[0])
                UT = int(row[1])

                lon1 = row[2]
                lat1 = row[3]
                h1 = row[4] / 1000

                lon2 = row[5]
                lat2 = row[6]
                h2 = row[7] / 1000

                time = NEQTime(mth, UT)
                NEQ_global = NequickG_global(time, BX)

                row[8] = NEQ_global.sTEC(h1, lat1, lon1, h2, lat2, lon2) / 10**16
                writer.writerow(row)


highsolaractivity()


# BX = GalileoBroadcast(121.129893, 0.351254133, 0.0134635348) # Medium solar activity
# BX = GalileoBroadcast(2.580271, 0.127628236, 0.0252748384) # Low solar activity


