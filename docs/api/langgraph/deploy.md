# LangSmith Deployment

[Production](/oss/python/langgraph/application-structure)

# LangSmith Deployment
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.This guide shows you how to deploy your agent to **[LangSmith Cloud](/langsmith/deploy-to-cloud)**, a fully managed hosting platform designed for agent workloads. With Cloud deployment, you can deploy directly from your GitHub repository—LangSmith handles the infrastructure, scaling, and operational concerns.
Traditional hosting platforms are built for stateless, short-lived web applications. LangSmith Cloud is **purpose-built for stateful, long-running agents** that require persistent state and background execution.
LangSmith offers multiple deployment options beyond Cloud, including deploying with a [control plane (hybrid/self-hosted)](/langsmith/deploy-with-control-plane) or as [standalone servers](/langsmith/deploy-standalone-server). For more information, refer to the [Deployment overview](/langsmith/deployment).

## [​](#prerequisites)Prerequisites

Before you begin, ensure you have the following:

{' ' * (self.list_depth - 1)}- A [GitHub account](https://github.com/)

{' ' * (self.list_depth - 1)}- A [LangSmith account](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langgraph-deploy) (free to sign up)

## [​](#deploy-your-agent)Deploy your agent

### [​](#1-create-a-repository-on-github)1. Create a repository on GitHub

Your application’s code must reside in a GitHub repository to be deployed on LangSmith. Both public and private repositories are supported. For this quickstart, first make sure your app is LangGraph-compatible by following the [local server setup guide](/oss/python/langgraph/studio#set-up-local-agent-server). Then, push your code to the repository.

### [​](#2-deploy-to-langsmith)2. Deploy to LangSmith

1[](#)

Navigate to LangSmith DeploymentLog in to [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langgraph-deploy). In the left sidebar, select **Deployments**.2[](#)

Create new deploymentClick the **+ New Deployment** button. A pane will open where you can fill in the required fields.3[](#)

Link repositoryIf you are a first time user or adding a private repository that has not been previously connected, click the **Add new account** button and follow the instructions to connect your GitHub account.4[](#)

Deploy repositorySelect your application’s repository. Click **Submit** to deploy. This may take about 15 minutes to complete. You can check the status in the **Deployment details** view.

### [​](#3-test-your-application-in-studio)3. Test your application in Studio

Once your application is deployed:

{' ' * (self.list_depth - 1)}- Select the deployment you just created to view more details.

{' ' * (self.list_depth - 1)}- Click the **Studio** button in the top right corner. Studio will open to display your graph.

### [​](#4-get-the-api-url-for-your-deployment)4. Get the API URL for your deployment

{' ' * (self.list_depth - 1)}- In the **Deployment details** view in LangGraph, click the **API URL** to copy it to your clipboard.

{' ' * (self.list_depth - 1)}- Click the `URL` to copy it to the clipboard.

### [​](#5-test-the-api)5. Test the API

You can now test the API:

{' ' * (self.list_depth - 1)}- Python
{' ' * (self.list_depth - 1)}- Rest API

{' ' * (self.list_depth - 1)}- Install LangGraph SDK:

```
pip install langgraph-sdk

```

{' ' * (self.list_depth - 1)}- Send a message to the agent:

```
from langgraph_sdk import get_sync_client # or get_client for async

client = get_sync_client(url="your-deployment-url", api_key="your-langsmith-api-key")

for chunk in client.runs.stream(
 None, # Threadless run
 "agent", # Name of agent. Defined in langgraph.json.
 input={
 "messages": [{
 "role": "human",
 "content": "What is LangGraph?",
 }],
 },
 stream_mode="updates",
):
 print(f"Receiving new event of type: {chunk.event}...")
 print(chunk.data)
 print("\n\n")

```

```
curl -s --request POST \
 --url <DEPLOYMENT_URL>/runs/stream \
 --header 'Content-Type: application/json' \
 --header "X-Api-Key: <LANGSMITH API KEY> \
 --data "{
 \"assistant_id\": \"agent\", `# Name of agent. Defined in langgraph.json.`
 \"input\": {
 \"messages\": [
 {
 \"role\": \"human\",
 \"content\": \"What is LangGraph?\"
 }
 ]
 },
 \"stream_mode\": \"updates\"
 }"

```

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/deploy.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Agent Chat UIPrevious](/oss/python/langgraph/ui)[LangSmith ObservabilityNext](/oss/python/langgraph/observability)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
