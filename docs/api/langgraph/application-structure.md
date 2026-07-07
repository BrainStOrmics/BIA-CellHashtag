# Application structure

[Production](/oss/python/langgraph/application-structure)

# Application structure
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.A LangGraph application consists of one or more graphs, a configuration file (`langgraph.json`), a file that specifies dependencies, and an optional `.env` file that specifies environment variables.
This guide shows a typical structure of an application and shows you how to provide the required configuration to deploy an application with [LangSmith Deployment](/langsmith/deployment).
LangSmith Deployment is a managed hosting platform for deploying and scaling LangGraph agents. It handles the infrastructure, scaling, and operational concerns so you can deploy your stateful, long-running agents directly from your repository. Learn more in the [Deployment documentation](/langsmith/deployment).

## [​](#key-concepts)Key concepts

To deploy using the LangSmith, the following information should be provided:

{' ' * (self.list_depth - 1)}- A [LangGraph configuration file](#configuration-file-concepts) (`langgraph.json`) that specifies the dependencies, graphs, and environment variables to use for the application.

{' ' * (self.list_depth - 1)}- The [graphs](#graphs) that implement the logic of the application.

{' ' * (self.list_depth - 1)}- A file that specifies [dependencies](#dependencies) required to run the application.

{' ' * (self.list_depth - 1)}- [Environment variables](#environment-variables) that are required for the application to run.

## [​](#file-structure)File structure

Below are examples of directory structures for applications:

{' ' * (self.list_depth - 1)}- Python (requirements.txt)
{' ' * (self.list_depth - 1)}- Python (pyproject.toml)
```
my-app/
├── my_agent # all project code lies within here
│ ├── utils # utilities for your graph
│ │ ├── __init__.py
│ │ ├── tools.py # tools for your graph
│ │ ├── nodes.py # node functions for your graph
│ │ └── state.py # state definition of your graph
│ ├── __init__.py
│ └── agent.py # code for constructing your graph
├── .env # environment variables
├── requirements.txt # package dependencies
└── langgraph.json # configuration file for LangGraph

```

```
my-app/
├── my_agent # all project code lies within here
│ ├── utils # utilities for your graph
│ │ ├── __init__.py
│ │ ├── tools.py # tools for your graph
│ │ ├── nodes.py # node functions for your graph
│ │ └── state.py # state definition of your graph
│ ├── __init__.py
│ └── agent.py # code for constructing your graph
├── .env # environment variables
├── langgraph.json # configuration file for LangGraph
└── pyproject.toml # dependencies for your project

```

The directory structure of a LangGraph application can vary depending on the programming language and the package manager used.
[]()

## [​](#configuration-file)Configuration file

The `langgraph.json` file is a JSON file that specifies the dependencies, graphs, environment variables, and other settings required to deploy a LangGraph application.
See the [LangGraph configuration file reference](/langsmith/cli#configuration-file) for details on all supported keys in the JSON file.
The [LangGraph CLI](/langsmith/cli) defaults to using the configuration file `langgraph.json` in the current directory.

### [​](#examples)Examples

{' ' * (self.list_depth - 1)}- The dependencies involve a custom local package and the `langchain_openai` package.

{' ' * (self.list_depth - 1)}- A single graph will be loaded from the file `./your_package/your_file.py` with the variable `variable`.

{' ' * (self.list_depth - 1)}- The environment variables are loaded from the `.env` file.

```
{
 "dependencies": ["langchain_openai", "./your_package"],
 "graphs": {
 "my_agent": "./your_package/your_file.py:agent"
 },
 "env": "./.env"
}

```

## [​](#dependencies)Dependencies

A LangGraph application may depend on other Python packages.
You will generally need to specify the following information for dependencies to be set up correctly:

{' ' * (self.list_depth - 1)}- 
A file in the directory that specifies the dependencies (e.g. `requirements.txt`, `pyproject.toml`, or `package.json`).

{' ' * (self.list_depth - 1)}- 
A `dependencies` key in the [LangGraph configuration file](#configuration-file-concepts) that specifies the dependencies required to run the LangGraph application.

{' ' * (self.list_depth - 1)}- 
Any additional binaries or system libraries can be specified using `dockerfile_lines` key in the [LangGraph configuration file](#configuration-file-concepts).

## [​](#graphs)Graphs

Use the `graphs` key in the [LangGraph configuration file](#configuration-file-concepts) to specify which graphs will be available in the deployed LangGraph application.
You can specify one or more graphs in the configuration file. Each graph is identified by a name (which should be unique) and a path for either: (1) the compiled graph or (2) a function that makes a graph is defined.

## [​](#environment-variables)Environment variables

If you’re working with a deployed LangGraph application locally, you can configure environment variables in the `env` key of the [LangGraph configuration file](#configuration-file-concepts).
For a production deployment, you will typically want to configure the environment variables in the deployment environment.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/application-structure.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[SubgraphsPrevious](/oss/python/langgraph/use-subgraphs)[TestNext](/oss/python/langgraph/test)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
