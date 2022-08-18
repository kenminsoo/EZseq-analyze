import matplotlib.pyplot as plt

string = ('protein binding cytoplasm biological regulation regulation of cellular process regulation of biological process cellular macromolecule metabolic process organonitrogen compound metabolic process cellular protein metabolic process organic substance biosynthetic process protein metabolic process biosynthetic process cytosol cellular biosynthetic process regulation of cellular metabolic process regulation of primary metabolic process macromolecule modification regulation of nitrogen compound metabolic process protein modification process cellular protein modification process positive regulation of biological process regulation of cellular process biological regulation regulation of biological process organonitrogen compound metabolic process cellular macromolecule metabolic process biosynthetic process protein metabolic process cellular biosynthetic process organic substance biosynthetic process regulation of cellular metabolic process cellular protein metabolic process response to stimulus regulation of nitrogen compound metabolic process regulation of primary metabolic process cellular response to stimulus positive regulation of biological process positive regulation of cellular process developmental process regulation of metabolic process macromolecule modification macromolecule biosynthetic process organelle organization protein modification process cellular protein modification process cellular nitrogen compound biosynthetic process')

split_string = string.split(" ")
skip_words = ['process', 'to', 'of', 'cellular', 'biological']

return_storage = {}
num_words = 0
for word in split_string:
    num_words += 1
    print(num_words)
    if word.lower() in skip_words:
        num_words = num_words - 1
        print(num_words)
    elif word.lower() in return_storage:
        return_storage[word.lower()] += 1
    else:
        return_storage[word.lower()] = 1

filtered_dict = dict(filter(lambda word_ct: (word_ct[1] > 5), return_storage.items()))

density_dict = {}
for key in filtered_dict:
    density_dict[key] = filtered_dict[key] / num_words

words = list(density_dict.keys())
values = density_dict.values()

plt.bar(range(len(density_dict)), values, tick_label = words)
plt.title('Density of Key Words in Ontology Analysis \nfor Breast Cancer Cell Lines')
plt.xlabel('Term')
plt.ylabel('Density of Occurances')
plt.show()

