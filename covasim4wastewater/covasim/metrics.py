import os
import dendropy


class Metrics:

    @staticmethod
    def _load_tree(path):
        """Load a tree from a file path."""
        return dendropy.Tree.get(
            path=path,
            schema="newick",
            preserve_underscores=True,
        )

    @staticmethod
    def _normalize_labels(tree):
        """
        Normalize tip labels so ground-truth trees and IQ-TREE trees match.
        Ground-truth labels look like: agent_1952|day_12|event_803
        IQ-TREE truncates at the first | producing: agent_1952
        This strips everything from the first | onward, in-place.
        """
        for leaf in tree.leaf_node_iter():
            leaf.taxon.label = leaf.taxon.label.split("|")[0]
        return tree

    @staticmethod
    def _put_on_shared_namespace(tree1, tree2):
        """
        RF and Faith's PD both require both trees to share the same
        TaxonNamespace. Re-parses both onto a fresh shared namespace,
        then normalizes labels on both so they can be compared.
        """
        tns = dendropy.TaxonNamespace()
        t1 = dendropy.Tree.get(
            data=tree1.as_string(schema="newick"),
            schema="newick",
            taxon_namespace=tns,
            preserve_underscores=True,
        )
        t2 = dendropy.Tree.get(
            data=tree2.as_string(schema="newick"),
            schema="newick",
            taxon_namespace=tns,
            preserve_underscores=True,
        )
        Metrics._normalize_labels(t1)
        Metrics._normalize_labels(t2)
        return t1, t2

    @staticmethod
    def count_mutations(tree_path):
        """
        Sum branch lengths across an entire tree.

        Parameters
        ----------
        tree_path : str   path to a Newick tree file (.nwk, .treefile, etc.)

        NOTE on units: gives an exact mutation count only if branch lengths
        are raw mutation counts (e.g. from TransmissionTree.construct_newick(
        branch_length_mode="count")). IQ-TREE branch lengths are in
        expected substitutions/site -- multiply by alignment length to
        approximate a mutation event count.

        Returns
        -------
        total : float
        """
        tree = Metrics._load_tree(tree_path)
        return sum(edge.length or 0 for edge in tree.preorder_edge_iter())

    @staticmethod
    def rf_distance(tree1_path, tree2_path, restrict_to_common_taxa=True):
        """
        Robinson-Foulds distance between two trees, pruned to shared taxa.

        Parameters
        ----------
        tree1_path, tree2_path : str   paths to Newick tree files

        Returns
        -------
        dict with rf, max_rf, normalized_rf, n_shared_taxa
        """
        t1, t2 = Metrics._put_on_shared_namespace(
            Metrics._load_tree(tree1_path),
            Metrics._load_tree(tree2_path),
        )

        if restrict_to_common_taxa:
            labels1 = {leaf.taxon.label for leaf in t1.leaf_node_iter()}
            labels2 = {leaf.taxon.label for leaf in t2.leaf_node_iter()}
            common  = labels1 & labels2

            if not common:
                raise ValueError("The two trees share no taxa in common.")

            if labels1 - common:
                t1.retain_taxa_with_labels(common)
            if labels2 - common:
                t2.retain_taxa_with_labels(common)

            n_taxa = len(common)
        else:
            n_taxa = len({leaf.taxon.label for leaf in t1.leaf_node_iter()})

        t1.encode_bipartitions()
        t2.encode_bipartitions()

        rf     = dendropy.calculate.treecompare.symmetric_difference(t1, t2)
        max_rf = 2 * (n_taxa - 3) if n_taxa > 3 else 1

        return {
            "rf":            rf,
            "max_rf":        max_rf,
            "normalized_rf": rf / max_rf if max_rf else 0.0,
            "n_shared_taxa": n_taxa,
        }

    @staticmethod
    def _scale_branch_lengths(tree, factor):
        """Multiply every branch length in the tree by factor, in-place."""
        for edge in tree.preorder_edge_iter():
            if edge.length is not None:
                edge.length *= factor
        return tree

    @staticmethod
    def _faith_pd(tree1_path, tree2_path, alignment_length=None):
        """
        Faith's Phylogenetic Diversity on the shared taxon set.
        Computed separately on each tree (restricted to shared tips),
        returning both PD values and their absolute difference.

        Parameters
        ----------
        tree1_path, tree2_path : str   paths to Newick tree files
        alignment_length : int or None
            If provided, IQ-TREE branch lengths (expected substitutions/site)
            on tree2 are multiplied by alignment_length to convert them to
            approximate mutation counts -- putting both trees on the same
            scale as the ground-truth tree (tree1), whose branch lengths are
            already raw mutation counts.
            If None, no scaling is applied (use only when both trees already
            share the same branch length units).

        Returns
        -------
        dict with pd_tree1, pd_tree2, pd_difference, n_shared_taxa
        """
        t1, t2 = Metrics._put_on_shared_namespace(
            Metrics._load_tree(tree1_path),
            Metrics._load_tree(tree2_path),
        )

        if alignment_length is not None:
            Metrics._scale_branch_lengths(t2, alignment_length)

        labels1 = {leaf.taxon.label for leaf in t1.leaf_node_iter()}
        labels2 = {leaf.taxon.label for leaf in t2.leaf_node_iter()}
        common  = labels1 & labels2

        if not common:
            raise ValueError("The two trees share no taxa in common.")

        if labels1 - common:
            t1.retain_taxa_with_labels(common)
        if labels2 - common:
            t2.retain_taxa_with_labels(common)

        pd1 = t1.length()
        pd2 = t2.length()

        return {
            "pd_tree1":        pd1,
            "pd_tree2":        pd2,
            "pd_difference":   abs(pd1 - pd2),
            "n_shared_taxa":   len(common),
            "alignment_length": alignment_length,
        }

    @staticmethod
    def tree_distance(tree1_path, tree2_path, metric="rf", alignment_length=None):
        """
        Compare two trees using the named metric.

        Parameters
        ----------
        tree1_path, tree2_path : str   paths to Newick tree files
        metric : "rf" or "faith_pd"
        alignment_length : int or None
            Required when metric="faith_pd" and tree2 is an IQ-TREE-inferred
            tree (branch lengths in expected substitutions/site). Multiplies
            tree2's branch lengths by this value to convert them to mutation
            counts, matching tree1's ground-truth units.
            Ignored when metric="rf" (RF is topology-only, no branch lengths).

        Returns
        -------
        dict (shape depends on metric)
        """
        if metric == "rf":
            return Metrics.rf_distance(tree1_path, tree2_path)
        elif metric == "faith_pd":
            return Metrics._faith_pd(tree1_path, tree2_path, alignment_length=alignment_length)
        else:
            raise ValueError(f"Unknown metric '{metric}'. Choose 'rf' or 'faith_pd'.")