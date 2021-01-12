import time_advancer as adv

def advance_classes():
    config = {}
    config["advancer"] = {}
    config["pop"] = {}
    config["model"] = {}

    # Advancer parameters
    config["advancer"]["method"] = "Euler"
    config["advancer"]["dt"] = 0.001
    config["advancer"]["T"] = 10
    config["advancer"]["p"] = 3

    # Population properties
    config["pop"]["NR0"] = 1
    config["pop"]["shape"] = "lognormal"
    config["pop"]["binning"] = "Simpson"
    config["pop"]["sigR0"] = 0.3

    # Bubble properties
    config["model"]["model"] = "RPE"
    config["model"]["R"] = 1.
    config["model"]["V"] = 0.
    config["model"]["gamma"] = 1.4
    # config["model"]["Ca"] = 0.5
    config["model"]["Re_inv"] = 0.001
    # config["model"]["Web"] = 10.

    myadv = adv.time_advancer(
            config=config["advancer"])
    myadv.initialize_state(
            pop_config=config["pop"],
            model_config=config["model"])
    myadv.run()

    # print(myadv.state.vals)
    # print(myadv.state.rhs)

    return

if __name__ == "__main__":

    advance_classes()
