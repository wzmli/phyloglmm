Dear Authors,

Thank you for submitting your research to MEE and apologies for the long time getting back to you, which partly a consequence of the ongoing lock-down work-from-home situation. Your manuscript has been reviewed by two highly qualified reviewers who found your work interesting while raising several points about the presentation and extent of the simulations, the scalability of the method and suggesting the inclusion of empirical examples. Based on the reviewers’ judgment and my own take on the paper I cannot recommend the acceptance of this paper at this point, but would encourage you to resubmit a thoroughly re-worked version.

MLi: Thank you for letting us resubmit.

As currently presented I found this paper a little bit in between an application note and a research paper but lacking some key components from both, but turning it into a research paper would require, in my opinion, the addition of an empirical example (see also comments from reviewer 1).

MLi: This paper is intended as an application/software paper. 
TODO: didn't we already do this? double check. We have included the example provided in Chapter 11 Garamszegi and Li et al. Dune mellow example. 

The benchmarking of the new method is interesting but is mixing too many confounding factors linked with the different implementations of the various packages (see also comment from reviewer 2). What is the performance of the new model all else being equal? How much of the performance results is actually due to differences the optimization among the different R packages?

MLi: The model is the mathematically equivalent, thus, the difference is the speed/efficiency.
MLi: Comparing pez and phyr uses lme4 correlation structure. This is the most direct comparision. 

In the Introduction, a broader perspective about PCMs should be given before jumping into the specifics of the models tested here. Phylogenetic comparative methods include lots of other things, see e.g. Luke Harmon’s book: https://lukejharmon.github.io/pcm/chapters/. In general, a better reference to previous work is required and some claims of the paper should be either demonstrated or removed. For example, stating that “existing [PCM] procedures are either insufficiently flexible or too computationally demanding to analyze large data sets” is simply not true. There are many research papers applying PCMs to large datasets (and as a side note: how large should a dataset be to be considered large?).

Todo: Revisit and rewrite. It is true many research papers used PCMs to large datasets, but it depends on the complexity of the model. I.e. If we simply use GLM methods in Table 1, it is relative fast even with large datasets. If we have more complicated models, then it will be very computationally demanding. 

A sentence like “Although many studies include multiple observations per species, phylogenetic analyses rarely take advantage of such information to partition variability more finely” needs references supporting it. For example, Kostikova et al. 2016 (DOI:10.1093/sysbio/syw010) present a model that seems relevant here. Further, claiming that ‘brms’ is more efficient than ‘MCMCglmm’ needs to be backed up by either a reference or a test (see also comments from the reviewers).

Todo: add reference!

In the methods it should be better clarified what is the state of the art and what are the novel parts you are presenting in this study.

Todo:... MCMCglmm?

In lines 74-75 “phylogenetic mixed models in community ecology” are mentioned but it is unclear what the link is to the previous paragraph and the to the following sentences.

Todo: revisit

Why are the parameter estimates of the two implementations of the new model different in Fig. 5? And, in the same figure, why were runs of alternative implementations killed after 30 min?

Todo: Help? lme4 vs glmmTMB differences? 

I am not sure about the point raised in the Discussion: “Establishing the practical level of model complexity for a given problem and data set is an open and difficult general problem throughout statistical modeling, not just in phylogenetic studies.” Model complexity is typically addressed by model testing, which is not mentioned here.

Todo: What do they mean by model testing here? This feels like snooping and cheating if we test and select models. 

In the discussion about within-species variance and measurement errors (l. 272-274) you could make several references to previous work, including several papers published in MEE.

Todo: find references.

I hope you will find these and the reviewers’ comments useful.
Best wishes,

Daniele Silvestro

Reviewer(s)' Comments to Author:
Reviewer: 1

Comments to the Corresponding Author
The authors propose a reformulation of the phylogenetic mixed model in order to make that class of PCM faster and more flexible. Given, their results, their goal has been reached. I found the paper well written and informative and their method of interest for a large audience. I don’t have any remarks about the scientific aspect nor about the benchmarking of that application which are thorough and well detailed.

MLi: Thanks.

However, I think the article could fail to reach the extent of its audience by missing several key elements:

- Because that work isn’t a theoretical work but instead the description of an applied statistics method, adding a few examples, such as the one provided at l68-73, would help the reader understanding the value of those models. I think examples should particularly be aimed at the “Simulation” § l132-155, where they could shed light on the subtle differences between the parameters
- In the same spirit of applied statistics, a brief description (or table) of which functions to use in the libraries lme4 and glmmTMB to make these models work would make it easier to start. There should also be mention of the tutorial available in the repository
- It would also be beneficial to add an application of the method with biological data directly in the article. The more your article is associated with biological ideas, the more users you will reach in my opinion.
- Naming that class of models (even though it is still pglmm) could help the reader follow the results section and identify the method as well.

These remarks might sound cosmetic but I think they might help to reach a broader audience and help all kind of users getting started with those models.

I recommend the authors to think on those ideas to improve the visibility of their method and the journal to ask them to do so in order to increase the impact of the article.

In a purely scientific perspective, I have no objections to publish that article.


Todo: Thanks, these are wonderful and helpful suggestion. We have included a tutorial in the released version on zenodo/github. We have now included the link in the paper.... 


Minor comments:

L26. In the Abstract, a claim has been made about the benefit of the method against unbalanced observational designs, but we lack insight on that point throughout the article. I think that sort of claim should be tested and I expect the Bayesian methods to perform better with that kind of problems.

Todo: Table 1 shows all the possbilities (double check) with existing package. 

L10. Small grammatical error: “can allow to incorporate”

Reviewer: 2

Comments to the Corresponding Author
This paper proposes an alternative formulation for the phylogenetic generalized linear mixed model which can utilize existing R packages for linear mixed models (lme4 and glmmTMB). The authors illustrate the advantages of their proposed method through simulations. The paper is well written and easy to follow. My main concern is that the proposed method is not scalable. That is, the computational efficiency of the proposed method only relies on the highly optimized implementation of lme4 and glmmTMB and its algorithm is actually "slower" than existing algorithms in phylogenetics. The followings are my specific comments:

Todo: Why is it slower? Do we need to profile the pieces? 

1. What is the computational complexity of the proposed method? Does the alternative formulation help avoid matrix inversion such that the complexity is linear with respect to the number of species? Note that this is the case for many phylogenetic packages including phylolm and MCMCglmm. I am aware that it is hard or even impossible to achieve linear-complexity for all models. However, could linear-complexity be achieved by the proposed methods for some simple models such as the standard phylogenetic linear regression (equation 1)?

Todo: Exactly, we are avoiding the matrix inversion step. 

2. In the paper, the authors replace the package MCMCglmm by the package brms. Although the implementation of MCMCglmm is not as efficient as brms, the complexity of MCMCglmm is linear with respect to the number of species while brms should be cubic (since it doesn't have a special method to avoid matrix inversion). It may be a good idea to add a simulation with a series of trees from dozens to thousands of species for showing the scalability of these methods (including MCMCglmm).

Todo: Will look into it. 

3. The models in simulations are not sufficiently described. What are the formulae for the phylogenetic random intercept variance, phylogenetic random slope variance, and covariance between phylogenetic random intercept and slope, and so on? Do they depend on the phylogenetic tree? The authors provide a clear explanation for how their proposed method can be applied to the standard phylogenetic regression (equation 1) but do not include it in the simulation section. I suggest the authors add this model to the simulation part for illustrating the applicability of their proposed method for this simple case.

Todo: what? go back to the paper and see if we say this clearly, this should be trival and clear. 

+++++
RESUBMISSION INSTRUCTIONS
Please note that resubmitting your manuscript does not guarantee eventual acceptance, and that your resubmission may be subject to re-review before a decision is rendered. Please also ensure that your altered manuscript still conforms to our word limit of 3000 for applications.

Once you have made the suggested changes, go to https://mc.manuscriptcentral.com/mee-besjournals and login to your Author Centre. Click on "Manuscripts with Decisions," and then click on "Create a Resubmission" located next to the manuscript number. Then, follow the steps for resubmitting your manuscript.

Because we are trying to facilitate timely publication of manuscripts submitted to Methods in Ecology and Evolution, your new manuscript should be uploaded within 12 weeks. The deadline for your resubmission is 05-Aug-2020. If it is not possible for you to submit your manuscript by that date, please get in touch with the editorial office, otherwise we will consider your paper as a completely new submission.
