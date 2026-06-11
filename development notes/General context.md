# Main context

This file contains an index of context files and code paths for AI agents.

## The project

Right now, AsymmeTree only allows to simulate two level evolution: species-gene reconciliation. In this project we will extend the framework to allow three level evolution: host-symbiont-gene reconciliation.

We will stick to the reconciliation simulation approach, but we add the host-symbiont reconciliation level. Furthermore, to model gene loss and gains, we will rely on game theory.

## Context files

- Theoretical plan of development:  [Math and Algorithms.md](Math and Algorithms.md).
- Skill for contributing in the asymmetree:   [skill contributing.md](skill contributing.md).
- Skill for implementation:  [skill implementation.md](skill implementation.md).
- Task already addressed:  [old tasks](old tasks).

## Implementation files

- The development directory is  [Example code/](Example code/), everything goes inside this folder. 
- Code paths and structure:   [code paths and structure.md](code paths and structure.md) 
- To-do list:  [implementation.md](implementation.md) 
- Bug solving:  [bug_reports.md](bug_reports.md) 

## Workflow rules

- AI notes and explanatory implementation notes belong in [implemented.md](implemented.md), not in [implementation.md](implementation.md).
- In [implementation.md](implementation.md), AI should only change checkbox states unless the user explicitly asks for another kind of edit.
- New tasks added to [implementation.md](implementation.md) should be ignored until the user explicitly asks to address them.

## Ignore

-  [prompt.md](prompt.md) 
