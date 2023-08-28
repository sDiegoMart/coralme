import pickle

def load_pickle_me_model(directory):
    with open(directory, "rb") as f:
        return pickle.load(f)

def save_pickle_me_model(me, directory):
    with open(directory, "wb") as f:
        pickle.dump(me, f)
