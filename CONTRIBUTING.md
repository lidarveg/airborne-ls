# Contributing to JRSRP Projects

First off, thank you for contributing! By participating in this project,
you help advance remote sensing science and ensure our tools are robust, transparent, and reproducible.

## Intellectual Property & Licensing

This project is led by the University of Queensland, and all contributions are governed by the **MIT License** (see the [LICENSE](LICENSE) file).

* **For  Government Employees:** When you contribute, the copyright of your specific changes remains with **The State in which they are employed**,
but you are licensing those changes to the project and the public under the MIT license.


## Getting Started

### 1. Set Up Your Identity
Please use your professional email address for all commits to ensure proper attribution back to your organization.

```bash
git config user.name "Your Name"
git config user.email "your.name@department.qld.gov.au"
```

## Pre-commit Hooks

Before pushing code, you must install our pre-commit hooks. These automatically check for accidentally included secrets, credentials, or large data files.

```bash
uv add --dev pre-commit
uv run pre-commit install
```

## AI-Assisted Coding
We encourage the use of AI tools (e.g., Claude, GitHub Copilot) to increase productivity,
provided the following rules are followed:

* No Sensitive Context: Never paste internal server URLs, credentials, or Protected PII into AI prompts.
* Human Review: You are responsible for the code you commit. All AI-generated logic must be reviewed and tested by you.
* Transparency: For significant AI-generated blocks, a brief note in the Merge Request description is appreciated.


## The Development Process

### Merge Requests (MRs)
Fork/Branch: Create a feature branch for your work.

Small Commits: Keep commits focused and descriptive.

The code is your original work or properly attributed.

No sensitive government data is included.

You are authorized by your manager to contribute to this open project.

### Code Style
We use Ruff for Python formatting. Please run `ruff check .` before submitting your MR. Our CI/CD pipeline will also check this automatically.
