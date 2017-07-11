import json
with open('data.txt', 'r') as outfile:
    j =json.load(outfile)
data = j['PKN']
for i in data:
    print(i['source'])