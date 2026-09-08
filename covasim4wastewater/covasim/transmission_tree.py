from collections import defaultdict, deque


INT_TO_NUC = {0: "A", 1: "C", 2: "G", 3: "T"}


class TransmissionTree:
    def __init__(self, tt, ref_seq):
        self.tt = tt
        self.ref_seq = ref_seq

        self.parent_of = {}        # event_id (tgt, date) -> (parent_event_id or None, branch_mutations)
        self.children_of = defaultdict(list)
        self.seeds = []             # root events, one per independent introduction
        self.genomes = {}           # event_id -> reconstructed sequence (list of chars)
        self.event_label = {}       # event_id -> "agent_{aid}|day_{date}|event_{i}"

        self.newick = None

    def build_tree(self):
        """Link transmission events into a tree and reconstruct every
        node's genome from the reference + accumulated branch mutations."""
        latest_event = {}

        # link events into parent/child relationships
        for t in self.tt.infection_log:
            tgt, src = t["target"], t["source"]
            date, muts = int(t["date"]), t["branch_mutations"]
            event_id = (tgt, date)

            if src is None:
                self.parent_of[event_id] = (None, muts)
                self.seeds.append(event_id)
            else:
                self.parent_of[event_id] = (latest_event.get(src), muts)

            latest_event[tgt] = event_id

        for event, (parent_event, _mut) in self.parent_of.items():
            if parent_event is not None:
                self.children_of[parent_event].append(event)

        # reconstruct root genomes
        for s in self.seeds:
            _, mut = self.parent_of[s]
            seq = self.ref_seq.copy()
            for pos, ref_idx, alt_idx in mut:
                assert seq[pos] == INT_TO_NUC[ref_idx]
                seq[pos] = INT_TO_NUC[alt_idx]
            self.genomes[s] = seq

        # BFS down from each root to reconstruct every descendant genome 
        seen = set(self.seeds)
        queue = deque(self.seeds)

        while queue:
            node = queue.popleft()
            for child in self.children_of[node]:
                if child in seen:
                    continue
                seen.add(child)

                seq = self.genomes[node].copy()
                _, mut = self.parent_of[child]
                for pos, ref_idx, alt_idx in mut:
                    seq[pos] = INT_TO_NUC[alt_idx]
                self.genomes[child] = seq
                queue.append(child)

        # assign labels so that fa export and the newick tree share the same label
        self.event_label = {
            event: f"agent_{event[0]}|day_{int(event[1])}|event_{i}"
            for i, event in enumerate(self.genomes)
        }

        return self

    def export_fasta(self, filepath):
        """Write every reconstructed genome to a FASTA file."""
        with open(filepath, "w") as f:
            for event, seq in self.genomes.items():
                f.write(f">{self.event_label[event]}\n{''.join(seq)}\n")
        print(f"Exported {len(self.genomes)} sequences to {filepath}")

    def construct_newick(self, branch_length_mode="count"):
        """
        Returns one Newick string per independent introduction (root).
        branch_length_mode: "count" -> branch length = number of true
                             mutations on that edge; "unit" -> all 1s.
        """
        def build(event):
            label = self.event_label[event]
            kids = self.children_of.get(event, [])
            if not kids:
                return label

            parts = []
            for kid in kids:
                _, mut = self.parent_of[kid]
                bl = len(mut) if branch_length_mode == "count" else 1
                parts.append(f"{build(kid)}:{bl}")
            return "(" + ",".join(parts) + ")" + label
        self.newick = [build(root) + ";" for root in self.seeds]

    def save_newick(self, save_path):
        if self.newick is not None and isinstance(self.newick, list):
            # combine multiple trees into a single tree 
            subtrees = [t.rstrip(";") for t in self.newick]
            combined_newick = "(" + ",".join(f"{t}:0" for t in subtrees) + ")Root;"
            with open(save_path, "w") as f:
                f.write(combined_newick)
            print(f"Saved combined tree to {save_path}")
        else:
            raise ValueError("Call tree construction before saving")


