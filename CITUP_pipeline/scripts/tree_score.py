import fileinput

for i, line in enumerate(fileinput.input()):
    if i == 3:
        fields = line.rstrip()
        print(fields)
        break