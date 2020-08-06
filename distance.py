import os
import csv
import math

class Distance(object):
    def __init__(self, file_name):
        self.debug = False
        self.file_name = file_name
        self.header = ['SC', 'NETWORK_ID', 'LAT', 'LON']
        self.neigh_prefix = 'NEIGH'
        self.data = []

    def get_output_filename(self):
        file_parts = os.path.splitext(self.file_name)
        return os.path.basename(file_parts[0]) + '.csv'

    def read_data(self):
        header = []

        if not os.path.exists(self.file_name):
            raise Exception("data file is missing: '{}'".format(self.file_name))

        with open(self.file_name, newline='') as f:
            # read header
            reader = csv.reader(f, delimiter=' ', skipinitialspace=True, quoting=csv.QUOTE_NONE)
            for row in reader:
                for name in row:
                    header.append(name.upper())
                break

            # check header
            if header != self.header:
                raise Exception("header line is missing or incorrect")

            # read data
            for row in reader:
                self.data.append({'SC': int(row[0]),
                                  'NETWORK_ID': row[1],
                                  'LAT': float(row[2]),
                                  'LON': float(row[3])})

        if self.debug:
            print('HEADER:', header)
            print('DATA:', self.data)

    def calculate_distance(self, location1, location2):
        lat1, lon1 = location1
        lat2, lon2 = location2
        radius = 6371 # radius of the earth

        dlat = math.radians(lat2 - lat1)
        dlon = math.radians(lon2 - lon1)

        #a = math.sin(dlat/2) * math.sin(dlat/2) + math.cos(math.radians(lat1)) \
        #    * math.cos(math.radians(lat2)) * math.sin(dlon/2) * math.sin(dlon/2)
        a = math.sin(dlat / 2) ** 2 \
            + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * (math.sin(dlon / 2) ** 2)
        c = 2 * math.atan2(math.sqrt(a), math.sqrt(1 - a))
        d = radius * c

        return d

    def get_distance_matrix_with_sc(self, sc):
        group_sc = list(filter(lambda x: sc == x['SC'], self.data))
        header_sc = []

        matrix = []
        for row in group_sc:
            matrix.append([None] * len(group_sc))
        i = 0
        j = 0

        for hor in group_sc:
            j = 0
            for ver in group_sc:
                if hor['NETWORK_ID'] == ver['NETWORK_ID']:
                    matrix[i][j] = 0.0
                    header_sc.append(hor['NETWORK_ID'])
                else:
                    if matrix[j][i] is None:
                        matrix[i][j] = self.calculate_distance((hor['LAT'], hor['LON']), (ver['LAT'], ver['LON']))
                    else: # is already calculated
                        matrix[i][j] = matrix[j][i]
                j = j + 1
            i = i + 1

        return (header_sc, matrix)

    # https://stackoverflow.com/a/38406223
    def format_matrix(self, header, matrix, top_format, left_format, cell_format, row_delim, col_delim):
        table = [[''] + header] + [[name] + row for name, row in zip(header, matrix)]
        table_format = [['{:^{}}'] + len(header) * [top_format]] \
                     + len(matrix) * [[left_format] + len(header) * [cell_format]]
        col_widths = [max(
                          len(format.format(cell, 0))
                          for format, cell in zip(col_format, col))
                      for col_format, col in zip(zip(*table_format), zip(*table))]
        return row_delim.join(
                   col_delim.join(
                       format.format(cell, width)
                       for format, cell, width in zip(row_format, row, col_widths))
                   for row_format, row in zip(table_format, table))

    def print_distance_with_sc(self, sc):
        header, matrix = self.get_distance_matrix_with_sc(sc)
        print(self.format_matrix(header,
                                 matrix,
                                 '{:^{}}', '{:<{}}', '{:>{}.3f}', '\n', ' | '))

    def save_data(self):
        i = 0
        max_neigh = 0
        matrix = []
        is_group_start = False
        for row in self.data:
            i = 0
            if not matrix:
                is_group_start = True
            if is_group_start:
                is_group_start = False
                sc = row['SC']
                _, matrix = self.get_distance_matrix_with_sc(sc)
                if self.debug:
                    print('MATRIX:', matrix)
                matrix = list(reversed(matrix))
            distances = matrix.pop()
            for distance in distances:
                i = i + 1
                if i > max_neigh:
                    max_neigh = i
                row[self.neigh_prefix.upper() + str(i)] = distance

        header = self.header
        for i in range(1, max_neigh + 1):
            header.append(self.neigh_prefix.upper() + str(i))

        if self.debug:
            print('NEW HEADER:', header)

        for row in self.data:
            diff = len(header) - len(row)
            if diff > 0:
                for i in range(diff - 1, -1, -1):
                    row[self.neigh_prefix + str(max_neigh - i)] = ''

        print(self.data)

        output = self.get_output_filename()
        if self.debug:
            print('OUTPUT FILE:', output)

        with open(output, 'w', newline='') as csvfile:
            fieldnames = header
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()

            for row in self.data:
                writer.writerow(row)

# execute only if run as a script
if __name__ == "__main__":
    distance = Distance('data.txt')
    distance.read_data()
    #distance.print_distance_with_sc(0)
    distance.save_data()


