import array
import numpy as np

def ctopy_read(name):
    file_txt = open(name+".txt", "r")
    file_bin = open(name+".bin", "rb")
    data_dict = {}
    for entry in file_txt:
        words = entry.split()
        var_name = words[0]
        var_type = words[1]
        num_dims = int(words[2])
        var_len = 1
        dims = []
        for k in range(num_dims):
            dim = int(words[3+k])
            dims.append(dim)
            var_len *= dim
        data = array.array(var_type)
        data.fromfile(file_bin, var_len)
        data = np.array(data)
        data = np.reshape(data, dims, order='C')
        data_dict[var_name] = data
    file_txt.close()
    file_bin.close()
    return data_dict

################################################################################
# Testing
################################################################################
if __name__=="__main__":
    imported_data = ctopy_read("my_data")
    keys = imported_data.keys()
    for key in keys:
        print("Variable Name: ", key, " has value: \n", imported_data[key])
        print("Its type is: ", type(imported_data[key]), "\n")
