# Panel selection (D-097)

Produced by `python3 analysis/panel_select.py` before any selected paper was read. Rule in the script's docstring and in D-097.

## Readouts, by circuits covered

| readout | filled now | covered | new papers | circuits covered |
|---|---|---|---|---|
| locomotor_activity **(panel)** | 3 | 11 | 40 | arc_agrp, arc_kisspeptin, arc_pomc, lha_mch, lha_orexin, lha_vgat, pvh_crh, pvh_oxt, scn_vip, tmn_histamine, vmh_sf1 |
| food_intake **(panel)** | 5 | 10 | 160 | arc_agrp, arc_kisspeptin, arc_pomc, lha_mch, lha_orexin, lha_vgat, pvh_crh, pvh_oxt, tmn_histamine, vmh_sf1 |
| wakefulness **(panel)** | 4 | 9 | 91 | arc_agrp, lha_mch, lha_orexin, lha_vgat, poa_galanin, pvh_crh, pvh_oxt, scn_vip, tmn_histamine |
| core_temperature **(panel)** | 1 | 9 | 20 | arc_agrp, arc_pomc, lha_orexin, lha_vgat, poa_galanin, pvh_crh, pvh_oxt, scn_vip, vmh_sf1 |
| body_weight **(panel)** | 3 | 7 | 70 | arc_agrp, arc_pomc, lha_mch, lha_orexin, lha_vgat, pvh_oxt, vmh_sf1 |
| energy_expenditure | 1 | 6 | 39 | arc_agrp, arc_pomc, lha_mch, lha_orexin, pvh_oxt, vmh_sf1 |
| corticosterone | 4 | 6 | 25 | arc_agrp, arc_pomc, lha_orexin, pvh_crh, pvh_oxt, scn_vip |
| blood_glucose | 2 | 6 | 23 | arc_agrp, arc_pomc, lha_orexin, pvh_crh, pvh_oxt, vmh_sf1 |
| heart_rate | 2 | 6 | 14 | arc_agrp, arc_pomc, lha_orexin, pvh_crh, pvh_oxt, scn_vip |

## The panel: filled now (F), new papers to read (number), nothing (-)

| circuit | locomotor_activity | food_intake | wakefulness | core_temperature | body_weight |
|---|---|---|---|---|---|
| arc_agrp | 3 | F | 1 | 3 | F |
| arc_kisspeptin | 1 | 2 | - | - | - |
| arc_pomc | 3 | F | - | F | 3 |
| lha_mch | 2 | F | 3 | - | F |
| lha_orexin | F | F | F | 2 | F |
| lha_vgat | 1 | 3 | 2 | 1 | 3 |
| poa_galanin | - | - | F | 1 | - |
| pvh_crh | F | 3 | F | 2 | - |
| pvh_oxt | 2 | 3 | 3 | 3 | 1 |
| scn_vip | F | - | 3 | 1 | - |
| tmn_histamine | 2 | 1 | F | - | - |
| vmh_sf1 | 1 | F | - | 2 | 3 |

60 cells: 16 filled now, 30 with papers to read, 14 with none.

## Reading set: 53 papers (at most 3 per empty cell)

* PMC8093934 (2021) -- CRISPR knockdown of Kcnq3 attenuates the M-current and increases excitability of NPY/AgRP neurons to alter ene -- for arc_agrp/locomotor_activity
* PMC10752887 (2023) -- Cannabis Sativa targets mediobasal hypothalamic neurons to stimulate appetite. -- for arc_agrp/locomotor_activity
* PMC9893607 (2023) -- Modulation of foraging-like behaviors by cholesterol-FGF19 axis. -- for arc_agrp/locomotor_activity
* PMC8089354 (2021) -- Chronic Stimulation of Arcuate Kiss1 Neurons Decreases Bone Mass in Female Mice -- for arc_kisspeptin/locomotor_activity
* PMC10005948 (2023) -- Increased intrinsic and synaptic excitability of hypothalamic POMC neurons underlies chronic stress-induced be -- for arc_pomc/locomotor_activity
* PMC10205441 (2023) -- TLR4 in POMC neurons regulates thermogenesis in a sex-dependent manner. -- for arc_pomc/locomotor_activity
* PMC9418835 (2022) -- <i>Magel2</i> knockdown in hypothalamic POMC neurons innervating the medial amygdala reduces susceptibility to -- for arc_pomc/locomotor_activity
* PMC12326619 (2025) -- Lack of glutamate neurotransmission in melanin-concentrating hormone neurons alters mouse reproduction and met -- for lha_mch/locomotor_activity
* PMC11259978 (2024) -- Glutamatergic signaling from melanin-concentrating hormone-producing neurons: A requirement for memory regulat -- for lha_mch/locomotor_activity, lha_mch/wakefulness
* PMC6617814 (2019) -- Effects of GABA and Leptin Receptor-Expressing Neurons in the Lateral Hypothalamus on Feeding, Locomotion, and -- for lha_vgat/locomotor_activity, lha_vgat/food_intake, lha_vgat/core_temperature, lha_vgat/body_weight
* PMC13722414 (2026) -- Oxytocin neurons signal state-dependent transitions from rest to thermogenesis and behavioral arousal. -- for pvh_oxt/locomotor_activity, pvh_oxt/wakefulness, pvh_oxt/core_temperature
* PMC11402739 (2024) -- Paraventricular oxytocin neurons impact energy intake and expenditure: projections to the bed nucleus of the s -- for pvh_oxt/locomotor_activity, pvh_oxt/food_intake, pvh_oxt/core_temperature
* PMC11564525 (2024) -- Histamine-tuned subicular circuit mediates alert-driven accelerated locomotion in mice. -- for tmn_histamine/locomotor_activity
* PMC7218123 (2020) -- Increased Sensitivity of Mice Lacking Extrasynaptic δ-Containing GABA<sub>A</sub> Receptors to Histamine Recep -- for tmn_histamine/locomotor_activity
* PMC4632173 (2015) -- Sex-dependent changes in metabolism and behavior, as well as reduced anxiety after eliminating ventromedial hy -- for vmh_sf1/locomotor_activity
* PMC13031403 (2026) -- Sex-specific differences in mediobasal hypothalamus in response to nutritional states. -- for arc_kisspeptin/food_intake
* PMC6103748 (2018) -- Estrogenic-dependent glutamatergic neurotransmission from kisspeptin neurons governs feeding circuits in femal -- for arc_kisspeptin/food_intake
* PMC11042726 (2024) -- Tesofensine, a novel antiobesity drug, silences GABAergic hypothalamic neurons. -- for lha_vgat/food_intake, lha_vgat/body_weight
* PMC7165479 (2020) -- Control of food approach and eating by a GABAergic projection from lateral hypothalamus to dorsal pons. -- for lha_vgat/food_intake, lha_vgat/body_weight
* PMC13677739 (2026) -- A hypothalamic circuit links high-fat diet exposure to anxiety-like behavior and anxiety-associated hyperphagi -- for pvh_crh/food_intake
* PMC13430618 (2026) -- Corticotropin-releasing factor neurons in the bed nucleus of the stria terminalis modulate avoidance behaviors -- for pvh_crh/food_intake
* PMC13045551 (2026) -- Norepinephrine neurons in the locus coeruleus and the nucleus of the solitary tract drive different stress-rel -- for pvh_crh/food_intake
* PMC13657252 (2026) -- Lateralized vagal oxytocin receptor signaling separately controls feeding and socioemotional functions via hyp -- for pvh_oxt/food_intake
* PMC12356352 (2025) -- Oxytocin neurons in the paraventricular and supraoptic hypothalamic nuclei bidirectionally modulate food intak -- for pvh_oxt/food_intake
* PMC9418084 (2022) -- Origin of thyrotropin-releasing hormone neurons that innervate the tuberomammillary nuclei. -- for tmn_histamine/food_intake
* PMC13732037 (2026) -- Orexin and pupillary modulation for hunger control. -- for arc_agrp/wakefulness
* PMC8275739 (2021) -- Lateral Hypothalamic Area Glutamatergic Neurons and Their Projections to the Lateral Habenula Modulate the Ane -- for lha_mch/wakefulness
* PMC9751402 (2022) -- Selective activation of the hypothalamic orexinergic but not melanin-concentrating hormone neurons following p -- for lha_mch/wakefulness
* PMC9459096 (2022) -- Whole-brain monosynaptic outputs and presynaptic inputs of GABAergic neurons in the vestibular nuclei complex  -- for lha_vgat/wakefulness
* PMC5617207 (2017) -- Neurochemical Heterogeneity Among Lateral Hypothalamic Hypocretin/Orexin and Melanin-Concentrating Hormone Neu -- for lha_vgat/wakefulness
* PMC12877408 (2026) -- Decreased Oxytocin Mediates PVN-CA2 and PVN-PrL in Sleep Deprivation-Induced Social Memory Deficits. -- for pvh_oxt/wakefulness
* PMC12376954 (2025) -- Oxytocin Receptor Expression and Activation in Parasympathetic Brainstem Cardiac Vagal Neurons. -- for pvh_oxt/wakefulness
* PMC12422183 (2025) -- Acute REM sleep deprivation alleviated depression-like behavior mediated by inhibiting VIP neurons in the mPFC -- for scn_vip/wakefulness
* PMC11905142 (2025) -- Oxytocin modulates inhibitory balance in the prelimbic cortex to support social memory consolidation during RE -- for scn_vip/wakefulness
* PMC12208536 (2025) -- Mitofusin 2 controls mitochondrial and synaptic dynamics of suprachiasmatic VIP neurons and related circadian  -- for scn_vip/wakefulness, scn_vip/core_temperature
* PMC10468093 (2023) -- GIRK2 potassium channels expressed by the AgRP neurons decrease adiposity and body weight in mice. -- for arc_agrp/core_temperature
* PMC10681932 (2023) -- High-fat diet feeding disrupts the coupling of thermoregulation to energy homeostasis. -- for arc_agrp/core_temperature
* PMC7837681 (2020) -- Cold-induced hyperphagia requires AgRP neuron activation in mice. -- for arc_agrp/core_temperature
* PMC9576062 (2022) -- The hypothalamic RFamide, QRFP, increases feeding and locomotor activity: The role of Gpr103 and orexin recept -- for lha_orexin/core_temperature
* PMC3589707 (2013) -- The physiological role of orexin/hypocretin neurons in the regulation of sleep/wakefulness and neuroendocrine  -- for lha_orexin/core_temperature
* PMC6868514 (2019) -- Galanin Neurons Unite Sleep Homeostasis and α2-Adrenergic Sedation. -- for poa_galanin/core_temperature
* PMC12279160 (2025) -- Amino Acid-Sensing Neurons in the Anterior Piriform Cortex Control Brown Adipose Tissue Thermogenesis. -- for pvh_crh/core_temperature
* PMC4911344 (2016) -- CRFR1 in AgRP Neurons Modulates Sympathetic Nervous System Activity to Adapt to Cold Stress and Fasting. -- for pvh_crh/core_temperature
* PMC12923983 (2026) -- Hypothalamic malate dehydrogenase 2 modulates systemic glucose metabolism through oxytocin-mediated thermogene -- for pvh_oxt/core_temperature
* PMC13229187 (2026) -- Projection-defined hypothalamic outputs differentially regulate thermogenesis and lipolysis. -- for vmh_sf1/core_temperature
* PMC11274921 (2024) -- Chemogenetic Excitation of Ventromedial Hypothalamic Steroidogenic Factor 1 (SF1) Neurons Increases Muscle The -- for vmh_sf1/core_temperature
* PMC13014027 (2026) -- Noradrenergic inhibition of definitive POMC neurons through direct and indirect mechanisms. -- for arc_pomc/body_weight
* PMC5323886 (2017) -- Cholinergic neurons in the dorsomedial hypothalamus regulate food intake. -- for arc_pomc/body_weight
* PMC7722761 (2020) -- Optogenetic stimulation of the liver-projecting melanocortinergic pathway promotes hepatic glucose production. -- for arc_pomc/body_weight
* PMC5490268 (2017) -- Caffeine inhibits hypothalamic A<sub>1</sub>R to excite oxytocin neuron and ameliorate dietary obesity in mice -- for pvh_oxt/body_weight
* PMC12076790 (2025) -- CPT1C deficiency in SF1 neurons impairs early metabolic adaptation to dietary fats, leading to obesity. -- for vmh_sf1/body_weight
* PMC11650306 (2025) -- TGR5 receptors in SF1-expressing neurons of the ventromedial hypothalamus regulate glucose homeostasis. -- for vmh_sf1/body_weight
* PMC8262364 (2021) -- Rap1 in the VMH regulates glucose homeostasis. -- for vmh_sf1/body_weight
