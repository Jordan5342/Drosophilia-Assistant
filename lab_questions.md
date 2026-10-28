# Lab questions for batch hypothesis runs

<!--
Edit freely, then run:  python run_questions.py
Only new or edited questions are rerun; unchanged ones keep their earlier results.

Format: one numbered question per line ("26. What converts ...").
Optional indented lines under a question override the automatic PubMed search:
    search: scribble Ras tumor invasion
    genes: scrib, Ras85D
Headings (## ...) and blank lines are ignored.
-->

## Neuro & behavior

1. How does internal state (hunger, mating status, sleep pressure) switch the fly between competing behaviors, and which neurons in the FlyWire connectome act as the switch?
    search: internal state behavioral choice neurons
2. Where does the whole-brain connectome fail to predict behavior? Which strongly connected pathways have no measurable effect when silenced?
    search: connectome behavior prediction neural circuits
3. What signal tells the brain that sleep need has been discharged, and does it come from neurons, glia or peripheral tissues?
    search: sleep homeostasis sleep need signal
4. How does one neuromodulator (dopamine or octopamine) reconfigure the same circuit to produce different outputs in different contexts?
    search: neuromodulation dopamine octopamine circuit
5. Why do isogenic flies raised in identical conditions show stable individual differences in behavior, and which developmental variation causes them?
    search: individuality behavior isogenic variation
6. What role do glia play in memory consolidation in the mushroom body, separate from their role in neuronal support?
    search: glia memory consolidation mushroom body
7. When light and temperature cues conflict, how do clock neurons decide which one sets the circadian phase?
    search: circadian clock neurons light temperature entrainment

## Development & growth

8. How does a growing organ tell the rest of the body that it has not yet reached its final size, and how is that signal turned off?
    search: organ size developmental timing Dilp8
    genes: Ilp8, Lgr3
9. How do morphogen gradients scale with tissue size so that pattern stays proportional in small and large flies?
    search: morphogen gradient scaling wing disc
    genes: dpp
10. What keeps developmental timing precise under nutritional or temperature stress, and which genes buffer it?
    search: developmental timing robustness ecdysone nutrition
11. How do mechanical forces during morphogenesis (e.g., germband extension, dorsal closure) feed back on gene expression?
    search: mechanical forces morphogenesis gene expression
12. Why do imaginal discs lose their capacity to regenerate late in larval development, and can that capacity be restored?
    search: imaginal disc regeneration capacity late larval
13. How do germline and intestinal stem cell niches adjust stem cell output to diet and tissue demand?
    search: stem cell niche diet intestinal germline

## Aging & metabolism

14. Which tissues are rate-limiting for lifespan, and can a tissue-specific intervention extend lifespan as much as a whole-body one?
    search: tissue-specific lifespan extension
15. Which components of dietary restriction (total protein, specific amino acids, sugar ratio) are necessary and sufficient for lifespan extension?
    search: dietary restriction amino acids lifespan
16. How does the gut report nutritional status to the brain to shape food choice, and which enteroendocrine signals carry it?
    search: enteroendocrine gut brain food choice
17. Is age-related intestinal barrier failure a cause of death or a marker of it, and what triggers it?
    search: intestinal barrier dysfunction aging death
18. What tissue-level mechanisms produce sex differences in lifespan and in the response to lifespan interventions?
    search: sex differences lifespan
19. How does the fat body integrate signals from several organs to control insulin-like peptide release from the brain?
    search: fat body insulin-like peptide secretion
    genes: Ilp2
20. How does reproduction trade off against somatic maintenance, and which signal from the germline sets that balance?
    search: reproduction lifespan trade-off germline

## Disease models

21. Which genetic modifiers suppress neurodegeneration across fly models of tau, alpha-synuclein and TDP-43, and which act in all three?
    search: neurodegeneration genetic modifiers tau TDP-43 synuclein
22. Can human variants of uncertain significance be classified at scale in humanized flies, and how well do fly results agree with clinical data?
    search: humanized fly variants of uncertain significance
23. Why are some neuron types selectively vulnerable in fly Parkinson's models while their neighbors survive?
    search: dopaminergic neuron selective vulnerability Parkinson
24. How does neuronal activity level interact with mitochondrial quality control (PINK1/parkin) to determine neuron survival?
    search: PINK1 parkin neuronal activity mitophagy
    genes: Pink1, park
25. Which tumor-secreted factors are necessary for wasting of distant tissues in fly cachexia models?
    search: tumor cachexia wasting secreted factors
26. What converts benign overgrowth into invasive growth when an oncogene is combined with loss of cell polarity?
    search: scribble Ras tumor invasion
    genes: scrib, Ras85D
27. Can flies predict which approved drugs rescue phenotypes of a specific patient variant, and how often do those predictions hold in patients?
    search: patient variant drug screening fly model

## Immunity & microbiome

28. How does the gut tell commensal bacteria from pathogens when both carry similar microbial patterns?
    search: gut immunity commensal pathogen discrimination
29. Which microbial metabolites mediate the microbiome's effects on fly development rate, lifespan and behavior?
    search: microbiota metabolites development lifespan
30. What is the cost of immune activation for reproduction and lifespan, and which tissue pays it?
    search: immune activation cost reproduction lifespan
31. Through what molecular mechanism does Wolbachia protect its host against RNA viruses?
    search: Wolbachia antiviral protection
32. Does the fly immune system show trained-immunity-like memory, and if so, in which cells is it stored?
    search: innate immune memory priming
33. How do antiviral RNAi and inducible immune pathways divide the work of controlling infection across tissues?
    search: antiviral RNAi immune pathways
    genes: Dcr-2

## Gene regulation & genome

34. Which features of enhancer sequence (motif grammar, spacing, affinity) best predict where and when a gene is expressed?
    search: enhancer grammar motif transcription
35. Can sequence-to-expression models trained on fly data correctly predict the in vivo effect of noncoding variants?
    search: sequence to expression model enhancer variants
36. How do young (de novo) genes become essential, and how many of them are required for viability or fertility?
    search: de novo genes essential fertility
37. What happens to transposable element control in somatic tissues with age, and does it contribute to aging?
    search: transposable elements aging somatic
38. Why does disrupting 3D genome structure often have modest effects on gene expression in flies, and when does it matter?
    search: 3D genome topology gene expression
39. How does the dosage compensation complex find the X chromosome so specifically, and what happens when targeting is partial?
    search: dosage compensation MSL complex X chromosome
    genes: msl-2

## Cell biology

40. What does cell competition actually measure as cell fitness, and how do loser cells sense they are losing?
    search: cell competition loser cell elimination
41. How is apical-basal polarity maintained during rapid proliferation in growing epithelia?
    search: apical basal polarity epithelial proliferation
    genes: scrib, dlg1
42. What decides whether a stressed cell dies or triggers compensatory proliferation in its neighbors?
    search: apoptosis induced compensatory proliferation
43. How is spindle orientation coupled to fate determinant segregation in neuroblasts, and what happens when the coupling is weakened?
    search: neuroblast asymmetric division spindle orientation
    genes: insc, mud
44. How do epithelial cells sense crowding and choose which cell is extruded?
    search: epithelial cell extrusion crowding
45. How do cells coordinate organelle positioning with signaling during polarized growth?
    search: organelle positioning polarized growth signaling

## Evolution & ecology

46. Which genetic changes underlie the rapid seasonal adaptation seen in wild Drosophila populations?
    search: seasonal adaptation wild populations
47. How did host specialization in D. sechellia evolve at the level of sensory neurons and central circuits?
    search: Drosophila sechellia host specialization olfactory
48. Why is a gene essential in one Drosophila species but dispensable in a close relative?
    search: gene essentiality species divergence
49. How does developmental temperature shape adult physiology and behavior, and which effects are adaptive?
    search: developmental temperature plasticity adult

## Methods & rigor

50. How often do classic RNAi or GAL4-based phenotypes reproduce when retested with clean CRISPR alleles in a matched background?
    search: RNAi off-target CRISPR mutant phenotype reproducibility
