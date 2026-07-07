Run the LATS MCTS algorithm to find the best cell type annotation.
1. Initialize root node with cluster markers
2. Iterate: select (UCB1), expand (generate hypotheses), validate (CellWiki), evaluate (score), backpropagate
3. Extract the best node as final annotation
Read the lats-search SKILL.md for protocol details.
Use cellwiki skill to query marker evidence.
Use core/mcts.py for MCTSNode, select_node, backpropagate, extract_best.