# Submitting climt to Conda-Forge

This guide explains how to submit `climt` to the `conda-forge` channel.

## Prerequisites

1.  A GitHub account.
2.  Basic knowledge of Git and GitHub.

## Steps

1.  **Fork staged-recipes**
    - Go to [conda-forge/staged-recipes](https://github.com/conda-forge/staged-recipes) and click "Fork" to create a copy in your account.

2.  **Clone your fork**
    ```bash
    git clone https://github.com/YOUR_USERNAME/staged-recipes.git
    cd staged-recipes
    ```

3.  **Create a new branch**
    ```bash
    git checkout -b add-climt
    ```

4.  **Add the recipe**
    - Create a directory for `climt` inside `recipes/`.
    ```bash
    mkdir recipes/climt
    ```
    - Copy the `recipe/meta.yaml` from this repository to `recipes/climt/meta.yaml` in your clone of staged-recipes.

5.  **Update the recipe**
    - The provided `recipe/meta.yaml` is configured for local builds (`path: ..`). For conda-forge, you must point to a stable release artifact (e.g., PyPI or GitHub release).
    - Edit `recipes/climt/meta.yaml`:
      - Remove the `path: ..` line.
      - Add the `url` and `sha256` of the source tarball.

      Example:
      ```yaml
      source:
        url: https://pypi.io/packages/source/c/climt/climt-{{ version }}.tar.gz
        sha256: <SHA256_HASH_OF_THE_TARBALL>
      ```
      - You can find the SHA256 hash by downloading the tarball from PyPI or using `openssl sha256 <file>`.

    - Update the `recipe-maintainers` list:
      - Replace `- your-github-username` with your actual GitHub username.

6.  **Commit and Push**
    ```bash
    git add recipes/climt/meta.yaml
    git commit -m "Add recipe for climt"
    git push origin add-climt
    ```

7.  **Open a Pull Request**
    - Go to your fork on GitHub and open a Pull Request against the `main` branch of `conda-forge/staged-recipes`.
    - Follow the checklist in the PR template.
    - The conda-forge bot will run checks. If there are linting errors, fix them.

8.  **Maintenance**
    - Once merged, a new repository `conda-forge/climt-feedstock` will be created.
    - Future updates to `climt` will be handled in that feedstock via pull requests (often automatically opened by the bot).

## Notes

- The `setup.py` has been modified to be robust for conda builds (handling compilers provided by environment variables).
- Ensure that the dependencies listed in `meta.yaml` match `setup.py` and `requirements.txt`.
