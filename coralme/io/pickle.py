import pickle

def load_pickle_me_model(path):
    with open(path, "rb") as infile:
        return pickle.load(infile)

def save_pickle_me_model(me, path):
    with open(path, "wb") as outfile:
        pickle.dump(me, outfile)
