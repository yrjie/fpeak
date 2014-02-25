import hmm_faster

hmm_model = hmm_faster.HMM()
hmm_model.set_states(['Healthy', 'Fever'])
#hmm_model.set_states(['Fever', 'Healthy'])
hmm_model.set_observations(['normal', 'cold', 'dizzy'])
hmm_model.randomize_matrices(seed = 19997)
Pi_matrix={'Healthy': 0.6, 'Fever': 0.4}
T_matrix={'Healthy' : {'Healthy': 0.7, 'Fever': 0.3},
	'Fever' : {'Healthy': 0.4, 'Fever': 0.6},
	}
E_matrix={
	'Healthy' : {'normal': 0.5, 'cold': 0.4, 'dizzy': 0.1},
	'Fever' : {'normal': 0.1, 'cold': 0.3, 'dizzy': 0.6},
	}
a=['normal', 'cold', 'dizzy']
hmm_model.set_initial_matrix(Pi_matrix)
hmm_model.set_transition_matrix(T_matrix)
hmm_model.set_emission_matrix(E_matrix)
a=['normal', 'cold', 'dizzy','normal', 'cold', 'dizzy','normal', 'cold', 'dizzy']
hmm_model.train(a, max_iteration=1000, delta=0.001)
print hmm_model.get_model()
result=hmm_model.viterbi(a)
print result
