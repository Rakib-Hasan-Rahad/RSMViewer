# RSMViewer: A PyMOL plugin for RNA structural motif visualization

Rakib Hasan Rahad¹, Smriti Pranjal¹, Nabila Shahnaz Khan², Shaojie Zhang², and Cuncong Zhong¹,\*

¹Department of Electrical Engineering and Computer Science, University of Kansas, Lawrence, KS, USA 66045
²Department of Computer Science, University of Central Florida, Orlando, FL, USA 32816
\*Corresponding author: cczhong@ku.edu

> This file reproduces the main manuscript text for reference alongside the code.
> The peer-reviewed version of record is published at *Bioinformatics*; please
> cite that version. A companion supplementary document that expands on the
> implementation is in [SUPPLEMENT.md](SUPPLEMENT.md).

---

## Abstract

RNA structural motifs (RSMs) are recurrent three-dimensional RNA structural
elements characterized by conserved three-dimensional geometries and often
associated with specific molecular functions. Their identification, validation,
and functional characterization frequently require labor-intensive manual
inspection. However, existing molecular visualization tools provide limited
support for motif-centric exploration, making large-scale visualization and
comparative analysis of RNA motifs inefficient and time-consuming. To address
this limitation, we present RSMViewer, an open-source PyMOL plugin that
integrates heterogeneous motif annotations and supports batch visualization and
efficient structural superimposition of RNA motif instances. RSMViewer
streamlines the exploration and comparative analysis of RNA structural motifs,
substantially reducing manual effort and facilitating downstream structural and
functional interpretation.

**Availability and Implementation:** RSMViewer is implemented in Python and is
freely available at https://doi.org/10.5281/zenodo.22097090. It is compatible
with PyMOL 2.x and later, on all major operating systems.

**Supplementary Information:** Supplementary data are available at *Bioinformatics*
online. See [SUPPLEMENT.md](SUPPLEMENT.md).

---

## Introduction

RNA structural motifs (RSMs) are recurrent structural elements characterized by
their conserved three-dimensional (3D) conformation (Leontis et al., 2006).
Their greater-than-expected recurrence, particularly across nonhomologous
regions of otherwise unrelated or evolutionarily distant RNAs, suggests that
these 3D conformations are subject to strong evolutionary constraints and often
play important structural or molecular roles (Bohdan et al., 2023; Hendrix et
al., 2005; Leontis et al., 2006). The functional importance of RNA structural
motifs is exemplified by the kink-turn motif (Klein et al., 2001) in U4atac
snRNA (Benoit-Pilven et al., 2020), which forms a specific recognition site for
the 15.5K protein and is essential for the assembly of the minor spliceosome.
Mutations disrupting conserved nucleotides within this motif substantially
impair 15.5K binding and minor-intron splicing and are associated with severe
developmental disorders, including Taybi–Linder syndrome (MOPD I), illustrating
how perturbation of a conserved RNA motif can directly lead to molecular
dysfunction and disease (Jafarifar et al., 2014).

Given their structural conservation and recurrent occurrence, RNA structural
motifs can naturally be organized into families whose members share
characteristic structural and functional features. A prominent example of
family-based RNA classification is Rfam (Griffiths-Jones et al., 2003), which
uses expert-curated sequence alignments and consensus secondary structures to
organize noncoding RNAs and structured RNA elements into families. Such
family-level organization facilitates comparative analysis by revealing the
range of sequence and secondary-structural variation tolerated within a family
while preserving its characteristic biological function. This information, in
turn, helps define family-specific inclusion and exclusion criteria for
subsequent large-scale computational searches (Cui et al., 2016).

The success of Rfam in organizing RNA families using conserved sequence and
secondary-structure information provides a useful paradigm for extending
family-based classification to RNA tertiary-structure motifs, or RNA 3D
structural motifs. Unlike secondary-structure motifs, RNA 3D motifs are
primarily classified by conserved spatial arrangements of nucleotides and their
interactions, which often provide a more direct connection to molecular function
(Batey et al., 1999; Hendrix et al., 2005; Moore, 1999). However, a centralized,
comprehensive, and extensively expert-curated resource for RNA 3D structural
motifs remains unavailable. Consequently, researchers often need to consult
multiple databases or annotation resources when analyzing RNA 3D motifs.
Existing RNA 3D motif resources have largely been developed using computational
approaches, and different resources frequently employ distinct representations
and criteria for motif identification. For example, RNA 3D motifs may be
characterized according to backbone trajectories (Apostolico et al., 2009),
overall geometric similarity (Sarver et al., 2008), or characteristic patterns
of canonical and noncanonical base interactions (Zhong et al., 2010). Different
resources may also adopt different classification schemes, ranging from
hierarchical organizations that incorporate both secondary- and
tertiary-structural characteristics (Petrov et al., 2013) to flatter
classification systems based primarily on 3D structural similarity (Djelloul and
Denise, 2008; Nabila et al., 2025; Zhong and Zhang, 2012). These differences can
result in partially overlapping, inconsistent, or conflicting annotations for
the same structural region. In practice, researchers may therefore need to
visually inspect individual motif instances to evaluate their three-dimensional
geometry, correct annotation errors, or reconcile discrepancies among different
resources (Baulin et al., 2025; Li and Chen, 2023; Miskiewicz et al., 2017;
Nabila et al., 2025; Zhang et al., 2021; Zurkowski et al., 2025), creating a
substantial burden for systematic and large-scale RNA motif analysis.

To address these challenges, we developed RSMViewer, a PyMOL (Schrödinger, 2015)
plugin for RNA structural motif visualization and analysis. RSMViewer enables
users to visualize thousands of motif instances using only a few commands. It
also consolidates heterogeneous motif annotation and classification sources and
allows users to flexibly select and visualize shared, source-specific, or
conflicting annotations through SQL-like commands. Finally, RSMViewer implements
a medoid-based superimposition heuristic to facilitate comparative visualization
of structural variation among related motif instances and to identify
structurally conserved cores that may be associated with molecular function.
Excluding the one-time download of online database resources, RSMViewer can
visualize thousands of motif instances or superimpose dozens of them within a
few seconds. Collectively, these features streamline RNA structural motif
visualization, annotation integration, and comparative analysis, thereby
facilitating systematic exploration of RNA structural motifs across
heterogeneous annotation resources.

## Motif Representation and Annotation Consolidation

RSMViewer represents each motif instance internally as a set of residues,
without imposing constraints on the number of RNA strands or molecular chains
involved. This flexible representation allows motifs defined according to
different criteria, including sequence, secondary structure, noncanonical
interaction patterns, and 3D geometry, as well as motifs spanning arbitrary
numbers of RNA strands or molecular chains, to be represented within a unified
framework. The residue sets obtained from external sources are subject to
redundancy filtering. Specifically, if the Jaccard index between two residue sets
is less than 60%, they are considered distinct structural fragments and are both
retained. Otherwise, they are treated as representing the same structural
fragment and are merged. This filtering approach retains distinct motif
instances sharing minor residue set overlap, which is frequently seen from
neighboring motif instances sharing closing canonical base pairs. After
redundant residue sets have been resolved, each remaining residue set is
assigned a unique internal identifier to facilitate subsequent motif selection
and visualization.

RSMViewer subsequently constructs a two-dimensional table indexed by the
filtered motif residue sets, with annotations from individual sources stored as
attributes of each residue set. When an annotation source employs a hierarchical
motif classification scheme, each level of the hierarchy is represented by a
dedicated column. For example, one column may represent the secondary-structure
context of a motif (e.g., internal loop, external loop, or stem), whereas
another may represent its tertiary-structure classification (e.g., kink-turn,
sarcin–ricin, or C-loop). Annotation sources employing a flat classification
scheme are represented by a single column. This design preserves the
organization and semantics of the original annotation sources, thereby
minimizing information loss and reducing potential biases introduced during
annotation consolidation. The resulting representation also enables simple and
efficient SQL-like querying, selection, and analysis of motif instances,
facilitating the identification of motifs consistently annotated across multiple
sources as well as those uniquely annotated by individual sources.

## Large-scale Motif Superimposition

RSMViewer further implements a heuristic approach for efficient superimposition
of multiple motif instances. Because PyMOL's built-in "super" and "align"
commands operate on two structural objects at a time, RSMViewer extends pairwise
alignment to multiple motif instances through a medoid-based strategy. Given a
set of motif instances, RSMViewer first performs pairwise structural
superimposition and records the resulting RMSDs. It then identifies the medoid
as the motif instance with the lowest average RMSD to all other instances. All
remaining motif instances are subsequently superimposed onto the medoid. This
strategy places all motif instances within a common reference frame while
maintaining practical computational efficiency.

> **Figure 1:** Examples of RSMViewer's applications. (A) Visualization of a
> single RNA structural motif instance. (B) Visualization of multiple instances
> within a motif family. (C) Visualization of all motif instances from multiple
> families of interest. Different colors represent different families. (D)
> Superimposition of multiple motif instances.

## RSMViewer Workflow

A typical RSMViewer analysis workflow comprises three main steps: (1) structure
loading, (2) annotation loading, and (3) motif visualization and analysis. The
first structure-loading step imports RNA structures from the Protein Data Bank
(PDB) (Berman, 2000), either through direct online retrieval or from
user-specified local directories (e.g., newly-resolved RNA structures). Both the
legacy PDB (.pdb) and current mmCIF (.cif) formats are supported. RSMViewer can
load either a single RNA structure or multiple structures simultaneously.

The second annotation-loading step imports motif annotations, including family
labels and constituent residues, from external resources such as curated
databases and outputs generated by model-based motif search tools. RSMViewer
currently supports automatic retrieval of motif instances from the RNA 3D Motif
Atlas (Petrov et al., 2013) and Rfam (Ontiveros-Palacios et al., 2025) through
their respective APIs. RSMViewer also supports outputs from the model-based
motif search tools FR3D (Sarver et al., 2008) and RNAMotifScanX (Zhong and
Zhang, 2015). For these tools, RSMViewer first attempts to retrieve precomputed
search results from their online repositories. If such results are unavailable,
RSMViewer automatically executes the standalone versions of FR3D and
RNAMotifScanX to generate the required annotations. Users can customize the
search parameters for both tools by editing a human-readable configuration file.
Retrieved online resources are cached locally after the initial download,
thereby avoiding repeated and potentially time-consuming online retrieval.

Finally, in the visualization and analysis step, RSMViewer enables users to
issue SQL-like queries to select and compare motif instances across different
annotation sources, thereby supporting a range of downstream analyses. RSMViewer
can visualize either individual motif instances, specified by their unique
internal identifiers, or sets of motif instances selected according to
user-defined criteria. Example RSMViewer commands for several typical downstream
analyses are shown below.

## Example Use Cases of RSMViewer

**Application 1: visualize individual motif instances (Figure 1A):**

```text
rmv_fetch 1S72                                       # loads RNA structure
rmv_db RNA3DMotifAtlas                               # select source and loads annotations
# selects RNA3DMotifAtlas-annotated SARCIN-RICIN (SR) motifs in 1S72
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas, as group_SR
rmv_list group_SR                                    # lists selected motifs (IDs and residue sets)
rmv_view [motif_id]                                  # visualizes the specified motif
```

**Application 2: visualize all instances within a motif family (Figure 1B):**

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas, as group_SR
rmv_view group_SR                                    # visualizes all motifs in the group
```

**Application 3: visualize all instances within multiple motif families of interest (Figure 1C):**

```text
rmv_fetch 1S72, 1FFK                                 # loads multiple RNA structures
rmv_db RNA3DMotifAtlas, Rfam                         # loads multiple annotation sources
# selects the commonly annotated SR motifs from 1S72 and 1FFK
rmv_select SARCIN-RICIN, 1S72 and 1FFK, RNA3DMotifAtlas and Rfam, as group_SR
# selects the commonly annotated KT, CL, and EL motifs from all loaded structures
rmv_select KT, all, RNA3DMotifAtlas and Rfam, as group_KT
rmv_select CL, all, RNA3DMotifAtlas and Rfam, as group_CL
rmv_select EL, all, RNA3DMotifAtlas and Rfam, as group_EL
rmv_view group_SR, group_KT, group_CL, group_EL
```

**Application 4: Superimposes all instances within a motif family to investigate within-class structural variation (Figure 1D):**

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas, as group_SR
rmv_create_object group_SR                           # creates PyMOL objects for all motifs in group_SR
rmv_super group_SR                                   # superimpose all motifs in the group
```

**Application 5: Benchmark the performance of a model-based search tool (e.g., RNAMotifScanX) when taking another database (e.g., RNA 3D Motif Atlas) as the ground truth reference.**

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas,RNAMotifScanX
# selects RNAMotifScanX-identified SR motifs with support (true positives)
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas and RNAMotifScanX, as group_TP
# selects RNAMotifScanX-identified SR motifs without support (false positives)
rmv_select SARCIN-RICIN, 1S72, not RNA3DMotifAtlas and RNAMotifScanX, as group_FP
# selects SR motifs potentially missed by RNAMotifScanX (false negatives)
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas and not RNAMotifScanX, as group_FN
```

**Application 6: Manually inspect false positives to identify potentially novel motif instances through structural superimposition.**

```text
rmv_fetch 1S72
rmv_db RNA3DMotifAtlas,RNAMotifScanX
# selects RNAMotifScanX-identified SR motifs without support (false positives)
rmv_select SARCIN-RICIN, 1S72, not RNA3DMotifAtlas and RNAMotifScanX, as group_FP
# selects known Sarcin-Ricin motif instances annotated in the RNA 3D Motif Atlas
rmv_select SARCIN-RICIN, 1S72, RNA3DMotifAtlas, as group_known
# assigns distinct colors to the two SR motif groups
rmv_set_color group_FP, red
rmv_set_color group_known, blue
# combines the two groups
rmv_combine group_FP, group_known, as group_combined
rmv_create_object group_combined
rmv_super group_combined
```

## Conclusion

In summary, RSMViewer provides a flexible and interactive platform for the
visualization and comparative analysis of RNA structural motifs through the
integration of heterogeneous annotation resources. RSMViewer consolidates
annotations from multiple sources while reducing redundancy and preserving
source-specific annotation diversity. Its medoid-based superimposition strategy
enables efficient structural comparison across multiple motif instances and
facilitates the identification of conserved structural cores and variable
geometric features. Collectively, RSMViewer bridges large-scale motif annotation
and structure-centric interpretation, providing a practical framework for
systematic investigation of RNA motif structure, variation, and potential
functional relevance.

## Funding

This work was supported by the National Institutes of Health under grant number
R01GM102515 and by the National Science Foundation under grant number
DBI-1943291.

## References

- Apostolico, A., et al. Finding 3D motifs in ribosomal RNA structures. *Nucleic Acids Res* 2009;37(4):e29.
- Batey, R.T., Rambo, R.P. and Doudna, J.A. Tertiary Motifs in RNA Structure and Folding. *Angew Chem Int Ed Engl* 1999;38(16):2326-2343.
- Baulin, E.F., et al. ARTEM: a method for RNA and DNA tertiary motif identification with backbone permutations. *Genome Biology* 2025;26(1).
- Benoit-Pilven, C., et al. Clinical interpretation of variants identified in RNU4ATAC, a non-coding spliceosomal gene. *PLoS One* 2020;15(7):e0235655.
- Berman, H.M. The Protein Data Bank. *Nucleic Acids Research* 2000;28(1):235-242.
- Bohdan, D.R., et al. A comprehensive survey of long-range tertiary interactions and motifs in non-coding RNA structures. *Nucleic Acids Research* 2023;51(16):8367-8382.
- Cui, X., et al. CMsearch: simultaneous exploration of protein sequence space and structure space improves not only protein homology detection but also protein structure prediction. *Bioinformatics* 2016;32(12):i332-i340.
- Djelloul, M. and Denise, A. Automated motif extraction and classification in RNA tertiary structures. *RNA* 2008;14(12):2489-2497.
- Griffiths-Jones, S., et al. Rfam: an RNA family database. *Nucleic Acids Res* 2003;31(1):439-441.
- Hendrix, D.K., Brenner, S.E. and Holbrook, S.R. RNA structural motifs: building blocks of a modular biomolecule. *Q Rev Biophys* 2005;38(3):221-243.
- Jafarifar, F., et al. Biochemical defects in minor spliceosome function in the developmental disorder MOPD I. *RNA* 2014;20(7):1078-1089.
- Klein, D.J., et al. The kink-turn: a new RNA secondary structure motif. *EMBO J* 2001;20(15):4214-4221.
- Leontis, N.B., Lescoute, A. and Westhof, E. The building blocks and motifs of RNA architecture. *Curr Opin Struct Biol* 2006;16(3):279-287.
- Li, J. and Chen, S.J. RNAJP: enhanced RNA 3D structure predictions with non-canonical interactions and global topology sampling. *Nucleic Acids Res* 2023;51(7):3341-3356.
- Miskiewicz, J., et al. Bioinformatics Study of Structural Patterns in Plant MicroRNA Precursors. *Biomed Res Int* 2017;2017:6783010.
- Moore, P.B. Structural motifs in RNA. *Annu Rev Biochem* 1999;68:287-300.
- Nabila, Md and Zhang, S. GINClus: RNA structural motif clustering using graph isomorphism network. *NAR Genomics and Bioinformatics* 2025;7(2).
- Ontiveros-Palacios, N., et al. Rfam 15: RNA families database in 2025. *Nucleic Acids Research* 2025;53(D1):D258-D267.
- Petrov, A.I., Zirbel, C.L. and Leontis, N.B. Automated classification of RNA 3D motifs and the RNA 3D Motif Atlas. *RNA* 2013;19(10):1327-1340.
- Sarver, M., et al. FR3D: finding local and composite recurrent structural motifs in RNA 3D structures. *J Math Biol* 2008;56(1-2):215-252.
- Schrödinger, L.L.C. The PyMOL Molecular Graphics System, Version 3.0. 2015.
- Zhang, D., Chen, S.J. and Zhou, R. Modeling Noncanonical RNA Base Pairs by a Coarse-Grained IsRNA2 Model. *J Phys Chem B* 2021;125(43):11907-11915.
- Zhong, C., Tang, H. and Zhang, S. RNAMotifScan: automatic identification of RNA structural motifs using secondary structural alignment. *Nucleic Acids Res* 2010;38(18):e176.
- Zhong, C. and Zhang, S. Clustering RNA structural motifs in ribosomal RNAs using secondary structural alignment. *Nucleic Acids Res* 2012;40(3):1307-1317.
- Zhong, C. and Zhang, S. RNAMotifScanX: a graph alignment approach for RNA structural motif identification. *RNA* 2015;21(3):333-346.
- Zurkowski, M., et al. Detecting polynucleotide motifs: Pentads, hexads, and beyond. *PLoS Comput Biol* 2025;21(10):e1013633.
