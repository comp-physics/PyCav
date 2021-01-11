import time_advancer as adv

def advance_classes():
    config = {}
    config["advancer"] = {}
    config["bub"] = {}

    # Advancer parameters
    config["advancer"]["method"] = "Euler"
    config["advancer"]["dt"] = 0.1
    config["advancer"]["T"] = 1.

    # Bubble properties
    config["bub"]["model"] = "RPE"
    config["bub"]["NR0"] = 1
    config["bub"]["shape"] = "lognormal"
    config["bub"]["binning"] = "Simpson"
    config["bub"]["sigR0"] = 0.3

    myadv = adv.time_advancer(config["advancer"])
    myadv.initialize_state(config["bub"])
    myadv.run()
    # print(myadv.state.vals)
    # print(myadv.state.rhs)

    return

if __name__ == "__main__":

    advance_classes()
