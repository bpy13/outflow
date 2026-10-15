---
name: stakeholder-persona
description: Who judges this repo and how. Used by the external-reviewer agent and when deciding how to explain or present the outflow package and app.
---

# Stakeholders

## The author (owner)
First author of Yu et al. 2019 (MNRAS). Knows the physics and the paper code in depth; is not
a professional programmer. Cares that:
- results reproduce the published Figs. 6–9 (differences invisible by eye);
- the code reads alongside the paper (equation numbers, paper symbols);
- the repo looks professional enough to share publicly.

## Readers from academia
Astrophysicists who read the paper and want to rerun or extend the model. Judge by: can I
install it in minutes, find where eq. N lives, change a parameter and trust the output? Expect
cgs units, citations, and honest statements of assumptions and limits.

## Readers from a programming background
Developers with no astrophysics. Judge by: clear structure, names, tests, CI, a README that
explains what the model does in plain words, an app they can click through without the paper.

## What would disappoint them
- Results that silently differ from the paper, or defaults that differ from it.
- Physics hidden in cryptic one-letter variables or one giant function.
- An app that freezes without saying it is computing, or fails with a stack trace.
