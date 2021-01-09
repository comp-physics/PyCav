import sys
import advancer as adv
import numpy as np

def advance_classes():
    config = {}
    config["bub"] = {}
    config["advancer"] = {}

    config["bub"]["governing_dynamics"] = "RPE"
    config["bub"]["num_R0_nodes"] = 5
    config["bub"]["method"] = "Simpson"

    config["advancer"]["method"] = "Euler"
    config["advancer"]["time_step"] = 1.0e-5
    config["advancer"]["final_time"] = 30.0
    # config["advancer"]["error_tol"] = 1.0e-5
    config["advancer"]["num_steps"] = 20000
    config["advancer"]["num_steps_print"] = 1000
    config["advancer"]["num_steps_write"] = 1000
    config["advancer"]["output_dir"] = "D/"
    config["advancer"]["output_id"] = "example_2D"
    config["advancer"]["write_to"] = "txt"

    advancer = adv.time_advancer(config)

    # Initial condition
    mu1 = 1.0
    sigma1 = 0.1
    advancer.initialize_state_gaussian_univar(mu1, sigma1)
    advancer.run()

    return

if __name__ == "__main__":

    np.set_printoptions(formatter={"float": "{: 0.4E}".format})
    advance_classes()
