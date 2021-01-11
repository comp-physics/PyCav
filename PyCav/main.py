import time_advancer as adv

def advance_classes():
    config = {}
    config["advancer"] = {}
    config["state"] = {}
    config["model"] = {}

    # Advancer parameters
    config["advancer"]["method"] = "Euler"
    config["advancer"]["dt"] = 0.1
    config["advancer"]["T"] = 1.

    # Bubble properties
    config["state"]["NR0"] = 1
    config["state"]["shape"] = "lognormal"
    config["state"]["binning"] = "Simpson"
    config["state"]["sigR0"] = 0.3

    config["model"]["model"] = "RPE"
    config["model"]["R"] = 1.
    config["model"]["V"] = 0.
    # config["model"]["gamma"] = 1.4
    # config["model"]["Ca"] = 1.
    # config["model"]["Re_inv"] = 0
    # config["model"]["Web"] = 0

    myadv = adv.time_advancer(config["advancer"])
    myadv.initialize_state(config["state"],config["model"])
    myadv.run()
    # print(myadv.state.vals)
    # print(myadv.state.rhs)

    return

if __name__ == "__main__":

    advance_classes()
