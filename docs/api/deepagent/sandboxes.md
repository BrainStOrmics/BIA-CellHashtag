# Sandboxes

[Core capabilities](/oss/python/deepagents/harness)

# Sandboxes
Copy page

Execute code in isolated environments with sandbox backendsCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Agents generate code, interact with filesystems, and run shell commands. Because we can’t predict what an agent might do, it’s important that its environment is isolated so it can’t access credentials, files, or the network. Sandboxes provide this isolation by creating a boundary between the agent’s execution environment and your host system.
In Deep Agents, **sandboxes are [backends](/oss/python/deepagents/backends)** that define the environment where the agent operates. Unlike other backends (State, Filesystem, Store) which only expose file operations, sandbox backends also give the agent an `execute` tool for running shell commands. When you configure a sandbox backend, the agent gets:

{' ' * (self.list_depth - 1)}- All standard filesystem tools (`ls`, `read_file`, `write_file`, `edit_file`, `glob`, `grep`)

{' ' * (self.list_depth - 1)}- The `execute` tool for running arbitrary shell commands in the sandbox

{' ' * (self.list_depth - 1)}- A secure boundary that protects your host system

## [​](#why-use-sandboxes)Why use sandboxes?

Sandboxes are used for security.
They let agents execute arbitrary code, access files, and use the network without compromising your credentials, local files, or host system.
This isolation is essential when agents run autonomously.
Sandboxes are especially useful for:

{' ' * (self.list_depth - 1)}- Coding agents: Agents that run autonomously can use shell, git, clone repositories (many providers offer native git APIs, e.g., [Daytona’s git operations](https://www.daytona.io/docs/en/git-operations/)), and run Docker-in-Docker for build and test pipelines

{' ' * (self.list_depth - 1)}- Data analysis agents—Load files, install data analysis libraries (pandas, numpy, etc.), run statistical calculations, and create outputs like PowerPoint presentations in a safe, isolated environment

**Using Deep Agents Code?** Deep Agents Code has built-in sandbox support via the `--sandbox` flag. See [Use remote sandboxes](/oss/python/deepagents/code/remote-sandboxes) for Deep Agents Code-specific setup, flags (`--sandbox-id`, `--sandbox-setup`), and examples.

## [​](#basic-usage)Basic usage

These examples assume you have already created a sandbox/devbox using the provider’s SDK and have credentials set up. For signup, authentication, and provider-specific lifecycle details, see [Available providers](#available-providers).

{' ' * (self.list_depth - 1)}- Modal
{' ' * (self.list_depth - 1)}- Runloop
{' ' * (self.list_depth - 1)}- Daytona
{' ' * (self.list_depth - 1)}- LangSmithpipuv
```
pip install langchain-modal

```

```
import modal
from deepagents import create_deep_agent
from langchain_anthropic import ChatAnthropic
from langchain_modal import ModalSandbox

app = modal.App.lookup("your-app")
modal_sandbox = modal.Sandbox.create(app=app)
backend = ModalSandbox(sandbox=modal_sandbox)

agent = create_deep_agent(
 model=ChatAnthropic(model="claude-sonnet-4-6"),
 system_prompt="You are a Python coding assistant with sandbox access.",
 backend=backend,
)
try:
 result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Create a small Python package and run pytest",
 }
 ]
 }
 )
finally:
 modal_sandbox.terminate()

```
pipuv
```
pip install langchain-runloop

```

```
import os

from deepagents import create_deep_agent
from langchain_anthropic import ChatAnthropic
from langchain_runloop import RunloopSandbox
from runloop_api_client import RunloopSDK

client = RunloopSDK(bearer_token=os.environ["RUNLOOP_API_KEY"])

devbox = client.devbox.create()
backend = RunloopSandbox(devbox=devbox)

agent = create_deep_agent(
 model=ChatAnthropic(model="claude-sonnet-4-6"),
 system_prompt="You are a Python coding assistant with sandbox access.",
 backend=backend,
)

try:
 result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Create a small Python package and run pytest",
 }
 ]
 }
 )
finally:
 devbox.shutdown()

```
pipuv
```
pip install langchain-daytona

```

```
from daytona import Daytona
from deepagents import create_deep_agent
from langchain_anthropic import ChatAnthropic
from langchain_daytona import DaytonaSandbox

sandbox = Daytona().create()
backend = DaytonaSandbox(sandbox=sandbox)

agent = create_deep_agent(
 model=ChatAnthropic(model="claude-sonnet-4-6"),
 system_prompt="You are a Python coding assistant with sandbox access.",
 backend=backend,
)

try:
 result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Create a small Python package and run pytest",
 }
 ]
 }
 )
finally:
 sandbox.stop()

```
LangSmith sandboxes are currently in private beta.pipuv
```
pip install "langsmith[sandbox]"

```

```
from deepagents import create_deep_agent
from deepagents.backends import LangSmithSandbox
from langchain_anthropic import ChatAnthropic
from langsmith.sandbox import SandboxClient

client = SandboxClient()
ls_sandbox = client.create_sandbox()
backend = LangSmithSandbox(sandbox=ls_sandbox)

agent = create_deep_agent(
 model=ChatAnthropic(model="claude-sonnet-4-6"),
 system_prompt="You are a Python coding assistant with sandbox access.",
 backend=backend,
)
try:
 result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Create a small Python package and run pytest",
 }
 ]
 }
 )
finally:
 client.delete_sandbox(ls_sandbox.name)

```

## [​](#available-providers)Available providers

For provider-specific setup, authentication, and lifecycle details, see [sandbox integrations](/oss/python/integrations/sandboxes).
Don’t see your provider? You can implement your own sandbox backend. See [Contributing a sandbox integration](/oss/python/contributing/integrations-langchain).

## [​](#lifecycle-and-scoping)Lifecycle and scoping

Most applications choose either one sandbox per [thread](/langsmith/use-threads) (thread-scoped) or one shared sandbox for every thread on the same [assistant](/langsmith/assistants) (assistant-scoped).
Sandboxes consume resources and cost money until they are shut down. Make sure you shut sandboxes down once they are no longer in use.
For the full lifecycle table, async [graph factory](/langsmith/graph-rebuild) notes, TTL behavior, LangGraph Deployment wiring, and client-side examples, see [Sandbox lifecycle](/oss/python/deepagents/going-to-production#lifecycle) in Going to production.

### [​](#thread-scoped-default)Thread-scoped (default)

Each conversation gets its own sandbox. The first run creates it; follow-up turns on the same thread reuse it. When the thread ends or the sandbox TTL expires, the environment goes away. Store the mapping with provider labels or metadata as in the following example so each run resolves to the same sandbox.
When users can return after idle time, configure a TTL on the sandbox so the provider deletes or archives idle environments automatically.

{' ' * (self.list_depth - 1)}- Python
{' ' * (self.list_depth - 1)}- TypeScriptagent.py
```
from daytona import CreateSandboxFromSnapshotParams, Daytona
from deepagents import create_deep_agent
from langchain_core.runnables import RunnableConfig
from langchain_daytona import DaytonaSandbox

client = Daytona()

async def agent(config: RunnableConfig):
 thread_id = config["configurable"]["thread_id"]
 try:
 sandbox = await client.find_one(labels={"thread_id": thread_id})
 except Exception:
 sandbox = await client.create(
 CreateSandboxFromSnapshotParams(
 labels={"thread_id": thread_id},
 auto_delete_interval=3600, # TTL: clean up when idle
 )
 )
 return create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=DaytonaSandbox(sandbox=sandbox)
 )

```
src/agent.ts
```
import { Daytona } from "@daytonaio/sdk";
import { DaytonaSandbox } from "@langchain/daytona";
import { createDeepAgent } from "deepagents";
import type { LangGraphRunnableConfig } from "@langchain/langgraph";

const client = new Daytona();

export async function agent(config: LangGraphRunnableConfig) {
 const threadId = config.configurable?.thread_id as string;
 let sandbox;
 try {
 sandbox = await client.findOne({ labels: { thread_id: threadId } });
 } catch {
 sandbox = await client.create({
 labels: { thread_id: threadId },
 autoDeleteInterval: 3600, // TTL: clean up when idle
 });
 }
 return createDeepAgent({
 model: "google_genai:gemini-3.5-flash",
 backend: await DaytonaSandbox.fromId(sandbox.id),
 });
}

```

### [​](#assistant-scoped)Assistant-scoped

Every thread on the same assistant reuses one sandbox. Files, installed packages, and cloned repositories persist across conversations.
Assistant-scoped sandboxes accumulate in-sandbox state over time. Configure a TTL with your sandbox provider, use snapshots to reset periodically, or implement cleanup logic so disk and memory do not grow without bound.

{' ' * (self.list_depth - 1)}- Python
{' ' * (self.list_depth - 1)}- TypeScriptagent.py
```
from daytona import CreateSandboxFromSnapshotParams, Daytona
from deepagents import create_deep_agent
from langchain_core.runnables import RunnableConfig
from langchain_daytona import DaytonaSandbox

client = Daytona()

async def agent(config: RunnableConfig):
 assistant_id = config["configurable"]["assistant_id"]
 try:
 sandbox = await client.find_one(labels={"assistant_id": assistant_id})
 except Exception:
 sandbox = await client.create(
 CreateSandboxFromSnapshotParams(labels={"assistant_id": assistant_id})
 )
 return create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=DaytonaSandbox(sandbox=sandbox)
 )

```
src/agent.ts
```
import { Daytona } from "@daytonaio/sdk";
import { DaytonaSandbox } from "@langchain/daytona";
import { createDeepAgent } from "deepagents";
import type { LangGraphRunnableConfig } from "@langchain/langgraph";

const client = new Daytona();

export async function agent(config: LangGraphRunnableConfig) {
 const assistantId = config.configurable?.assistant_id as string;
 let sandbox;
 try {
 sandbox = await client.findOne({ labels: { assistant_id: assistantId } });
 } catch {
 sandbox = await client.create({ labels: { assistant_id: assistantId } });
 }
 return createDeepAgent({
 model: "google_genai:gemini-3.5-flash",
 backend: await DaytonaSandbox.fromId(sandbox.id),
 });
}

```

For manual create, execute, and teardown outside a graph factory, see [Basic usage](#basic-usage) and [sandbox integrations](/oss/python/integrations/sandboxes) for provider-specific APIs.

## [​](#integration-patterns)Integration patterns

There are two architecture patterns for integrating agents with sandboxes, based on where the agent runs.

### [​](#agent-in-sandbox-pattern)Agent in sandbox pattern

The agent runs inside the sandbox and you communicate with it over the network. You build a Docker or VM image with your agent framework pre-installed, run it inside the sandbox, and connect from outside to send messages.
Benefits:

{' ' * (self.list_depth - 1)}- ✅ Mirrors local development closely.

{' ' * (self.list_depth - 1)}- ✅ Tight coupling between agent and environment.

Trade-offs:

{' ' * (self.list_depth - 1)}- 🔴 API keys must live inside the sandbox (security risk).

{' ' * (self.list_depth - 1)}- 🔴 Updates require rebuilding images.

{' ' * (self.list_depth - 1)}- 🔴 Requires infrastructure for communication (WebSocket or HTTP layer).

To run an agent in a sandbox, build an image and install deepagents on it.

```
FROM python:3.11
RUN pip install deepagents-code

```

Then run the agent inside the sandbox.
To use the agent inside the sandbox you have to add additional infrastructure to handle communication between your application and the agent inside the sandbox.

### [​](#sandbox-as-tool-pattern)Sandbox as tool pattern

The agent runs on your machine or server. When it needs to execute code, it calls sandbox tools (such as `execute`, `read_file`, or `write_file`) which invoke the provider’s APIs to run operations in a remote sandbox.
Benefits:

{' ' * (self.list_depth - 1)}- ✅ Update agent code instantly without rebuilding images.

{' ' * (self.list_depth - 1)}- ✅ Cleaner separation between agent state and execution.

{' ' * (self.list_depth - 1)}- API keys stay outside the sandbox.

{' ' * (self.list_depth - 1)}- Sandbox failures don’t lose agent state.

{' ' * (self.list_depth - 1)}- Option to run tasks in multiple sandboxes in parallel.

{' ' * (self.list_depth - 1)}- ✅ Pay only for execution time.

Trade-offs:

{' ' * (self.list_depth - 1)}- 🔴 Network latency on each execution call.

Example
```
from daytona import Daytona
from deepagents import create_deep_agent
from dotenv import load_dotenv
from langchain_daytona import DaytonaSandbox

load_dotenv()

# Can also do this with AgentCore, E2B, Runloop, Modal
sandbox = Daytona().create()
backend = DaytonaSandbox(sandbox=sandbox)

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 backend=backend,
 system_prompt="You are a coding assistant with sandbox access. You can create and run code in the sandbox.",
)

try:
 result = agent.invoke(
 {
 "messages": [
 {
 "role": "user",
 "content": "Create a hello world Python script and run it",
 }
 ]
 }
 )
 print(result["messages"][-1].content)
except Exception:
 # Optional: delete the sandbox proactively on an exception
 sandbox.stop()
 raise

```

The examples in this doc use the sandbox as a tool pattern.
Choose the agent in sandbox pattern when your provider’s SDK handles the communication layer and you want production to mirror local development.
Choose the sandbox as tool pattern when you need to iterate quickly on agent logic, keep API keys outside the sandbox, or prefer cleaner separation of concerns.

## [​](#how-sandboxes-work)How sandboxes work

### [​](#isolation-boundaries)Isolation boundaries

All sandbox providers protect your host system from the agent’s filesystem and shell operations. The agent cannot read your local files, access environment variables on your machine, or interfere with other processes. However, sandboxes alone do **not** protect against:

{' ' * (self.list_depth - 1)}- **Context injection**: An attacker who controls part of the agent’s input can instruct it to run arbitrary commands inside the sandbox. The sandbox is isolated, but the agent has full control within it.

{' ' * (self.list_depth - 1)}- **Network exfiltration**: Unless network access is blocked, a context-injected agent can send data out of the sandbox over HTTP or DNS. Some providers support blocking network access (e.g., `blockNetwork: true` on Modal).

See [security considerations](#security-considerations) for how to handle secrets and mitigate these risks.

### [​](#the-execute-method)The `execute` method

Sandbox backends have a simple architecture: the only method a provider must implement is `execute()`, which runs a shell command and returns its output. Every other filesystem operation (`read`, `write`, `edit`, `ls`, `glob`, `grep`) is built on top of `execute()` by the [`BaseSandbox`](https://reference.langchain.com/python/deepagents/backends/sandbox/BaseSandbox) base class, which constructs scripts and runs them inside the sandbox via `execute()`.

This design means:

{' ' * (self.list_depth - 1)}- **Adding a new provider is straightforward.** Implement `execute()`—the base class handles everything else.

{' ' * (self.list_depth - 1)}- **The `execute` tool is conditionally available.** On every model call, the harness checks whether the backend implements [`SandboxBackendProtocol`](https://reference.langchain.com/python/deepagents/backends/protocol/SandboxBackendProtocol). If not, the tool is filtered out and the agent never sees it.

When the agent calls the `execute` tool, it provides a `command` string and gets back the combined stdout/stderr, exit code, and a truncation notice if the output was too large.
You can also call the backend `execute()` method directly in your application code.

{' ' * (self.list_depth - 1)}- Daytona
{' ' * (self.list_depth - 1)}- Modal
{' ' * (self.list_depth - 1)}- Runloop
{' ' * (self.list_depth - 1)}- AgentCore
{' ' * (self.list_depth - 1)}- LangSmithpipuv
```
pip install langchain-daytona

```

```
from daytona import Daytona

from langchain_daytona import DaytonaSandbox

sandbox = Daytona().create()
backend = DaytonaSandbox(sandbox=sandbox)

result = backend.execute("python --version")
print(result.output)

```

```
import modal

from langchain_modal import ModalSandbox

app = modal.App.lookup("your-app")
modal_sandbox = modal.Sandbox.create(app=app)
backend = ModalSandbox(sandbox=modal_sandbox)

result = backend.execute("python --version")
print(result.output)

```
pipuv
```
pip install langchain-runloop

```

```
from runloop_api_client import RunloopSDK

from langchain_runloop import RunloopSandbox

api_key = "..."
client = RunloopSDK(bearer_token=api_key)

devbox = client.devbox.create()
backend = RunloopSandbox(devbox=devbox)

try:
 result = backend.execute("python --version")
 print(result.output)
finally:
 devbox.shutdown()

```
pipuv
```
pip install langchain-agentcore-codeinterpreter

```

```
from bedrock_agentcore.tools.code_interpreter_client import CodeInterpreter

from langchain_agentcore_codeinterpreter import AgentCoreSandbox

interpreter = CodeInterpreter(region="us-west-2")
interpreter.start()

backend = AgentCoreSandbox(interpreter=interpreter)

try:
 result = backend.execute("python3 --version")
 print(result.output)
finally:
 interpreter.stop()

```

```
from langsmith.sandbox import SandboxClient

from deepagents.backends.langsmith import LangSmithSandbox

client = SandboxClient()
ls_sandbox = client.create_sandbox(template_name="deepagents-deploy")
backend = LangSmithSandbox(sandbox=ls_sandbox)

result = backend.execute("python --version")
print(result.output)

```

For example:

```
4
[Command succeeded with exit code 0]

```

```
bash: foobar: command not found
[Command failed with exit code 127]

```

If a command produces very large output, the result is automatically saved to a file and the agent is instructed to use `read_file` to access it incrementally. This prevents context window overflow.

### [​](#two-planes-of-file-access)Two planes of file access

There are two distinct ways files move in and out of a sandbox, and it’s important to understand when to use each:
**Agent filesystem tools**: `read_file`, `write_file`, `edit_file`, `ls`, `glob`, `grep`, and `execute` are the tools the LLM calls during its execution. These go through `execute()` inside the sandbox. The agent uses them to read code, write files, and run commands as part of its task.
**File transfer APIs**: the `uploadFiles()` and `downloadFiles()` methods that your application code calls. These use the provider’s native file transfer APIs (not shell commands) and are designed for moving files between your host environment and the sandbox. Use these to:

{' ' * (self.list_depth - 1)}- **Seed the sandbox** with source code, configuration, or data before the agent runs

{' ' * (self.list_depth - 1)}- **Retrieve artifacts** (generated code, build outputs, reports) after the agent finishes

{' ' * (self.list_depth - 1)}- **Pre-populate dependencies** that the agent will need

## [​](#working-with-files)Working with files

The deepagents sandbox backends support file transfer APIs for moving files between your application and the sandbox.

### [​](#seeding-the-sandbox)Seeding the sandbox

Use `upload_files()` to populate the sandbox before the agent runs. Paths must be absolute and contents are `bytes`:

{' ' * (self.list_depth - 1)}- Daytona
{' ' * (self.list_depth - 1)}- Modal
{' ' * (self.list_depth - 1)}- Runloop
{' ' * (self.list_depth - 1)}- AgentCore
{' ' * (self.list_depth - 1)}- LangSmithpipuv
```
pip install langchain-daytona

```

```
from daytona import Daytona

from langchain_daytona import DaytonaSandbox

sandbox = Daytona().create()
backend = DaytonaSandbox(sandbox=sandbox)

backend.upload_files(
 [
 ("/src/index.py", b"print('Hello')\n"),
 ("/pyproject.toml", b"[project]\nname = 'my-app'\n"),
 ]
)

```

```
import modal

from langchain_modal import ModalSandbox

app = modal.App.lookup("your-app")
modal_sandbox = modal.Sandbox.create(app=app)
backend = ModalSandbox(sandbox=modal_sandbox)

backend.upload_files(
 [
 ("/src/index.py", b"print('Hello')\n"),
 ("/pyproject.toml", b"[project]\nname = 'my-app'\n"),
 ]
)

```
pipuv
```
pip install langchain-runloop

```

```
from runloop_api_client import RunloopSDK

from langchain_runloop import RunloopSandbox

api_key = "..."
client = RunloopSDK(bearer_token=api_key)

devbox = client.devbox.create()
backend = RunloopSandbox(devbox=devbox)

backend.upload_files(
 [
 ("/src/index.py", b"print('Hello')\n"),
 ("/pyproject.toml", b"[project]\nname = 'my-app'\n"),
 ]
)

```
pipuv
```
pip install langchain-agentcore-codeinterpreter

```

```
from bedrock_agentcore.tools.code_interpreter_client import CodeInterpreter

from langchain_agentcore_codeinterpreter import AgentCoreSandbox

interpreter = CodeInterpreter(region="us-west-2")
interpreter.start()

backend = AgentCoreSandbox(interpreter=interpreter)

backend.upload_files(
 [
 ("hello.py", b"print('Hello')\n"),
 ("data.csv", b"name,value\na,1\nb,2\n"),
 ]
)

```

```
from langsmith.sandbox import SandboxClient

from deepagents.backends.langsmith import LangSmithSandbox

client = SandboxClient()
ls_sandbox = client.create_sandbox(template_name="deepagents-deploy")
backend = LangSmithSandbox(sandbox=ls_sandbox)

backend.upload_files(
 [
 ("/src/index.py", b"print('Hello')\n"),
 ("/pyproject.toml", b"[project]\nname = 'my-app'\n"),
 ]
)

```

### [​](#retrieving-artifacts)Retrieving artifacts

Use `download_files()` to retrieve files from the sandbox after the agent finishes:

{' ' * (self.list_depth - 1)}- Daytona
{' ' * (self.list_depth - 1)}- Modal
{' ' * (self.list_depth - 1)}- Runloop
{' ' * (self.list_depth - 1)}- AgentCore
{' ' * (self.list_depth - 1)}- LangSmithpipuv
```
pip install langchain-daytona

```

```
from daytona import Daytona

from langchain_daytona import DaytonaSandbox

sandbox = Daytona().create()
backend = DaytonaSandbox(sandbox=sandbox)

results = backend.download_files(["/src/index.py", "/output.txt"])
for result in results:
 if result.content is not None:
 print(f"{result.path}: {result.content.decode()}")
 else:
 print(f"Failed to download {result.path}: {result.error}")

```

```
import modal

from langchain_modal import ModalSandbox

app = modal.App.lookup("your-app")
modal_sandbox = modal.Sandbox.create(app=app)
backend = ModalSandbox(sandbox=modal_sandbox)

results = backend.download_files(["/src/index.py", "/output.txt"])
for result in results:
 if result.content is not None:
 print(f"{result.path}: {result.content.decode()}")
 else:
 print(f"Failed to download {result.path}: {result.error}")

```
pipuv
```
pip install langchain-runloop

```

```
from runloop_api_client import RunloopSDK

from langchain_runloop import RunloopSandbox

api_key = "..."
client = RunloopSDK(bearer_token=api_key)

devbox = client.devbox.create()
backend = RunloopSandbox(devbox=devbox)

results = backend.download_files(["/src/index.py", "/output.txt"])
for result in results:
 if result.content is not None:
 print(f"{result.path}: {result.content.decode()}")
 else:
 print(f"Failed to download {result.path}: {result.error}")

```
pipuv
```
pip install langchain-agentcore-codeinterpreter

```

```
from bedrock_agentcore.tools.code_interpreter_client import CodeInterpreter

from langchain_agentcore_codeinterpreter import AgentCoreSandbox

interpreter = CodeInterpreter(region="us-west-2")
interpreter.start()

backend = AgentCoreSandbox(interpreter=interpreter)

results = backend.download_files(["hello.py"])
for result in results:
 if result.content is not None:
 print(f"{result.path}: {result.content.decode()}")
 else:
 print(f"Failed to download {result.path}: {result.error}")

interpreter.stop()

```

```
from langsmith.sandbox import SandboxClient

from deepagents.backends.langsmith import LangSmithSandbox

client = SandboxClient()
ls_sandbox = client.create_sandbox(template_name="deepagents-deploy")
backend = LangSmithSandbox(sandbox=ls_sandbox)

results = backend.download_files(["/src/index.py", "/output.txt"])
for result in results:
 if result.content is not None:
 print(f"{result.path}: {result.content.decode()}")
 else:
 print(f"Failed to download {result.path}: {result.error}")

```

Inside the sandbox, the agent uses filesystem tools (`read_file`, `write_file`). The `upload_files` and `download_files` methods are for your application code to move files across the boundary between your host and the sandbox.

## [​](#security-considerations)Security considerations

Sandboxes isolate code execution from your host system, but they don’t protect against **context injection**. An attacker who controls part of the agent’s input can instruct it to read files, run commands, or exfiltrate data from within the sandbox. This makes credentials inside the sandbox especially dangerous.
**Never put secrets inside a sandbox.** API keys, tokens, database credentials, and other secrets injected into a sandbox (via environment variables, mounted files, or the `secrets` option) can be read and exfiltrated by a context-injected agent. This applies even to short-lived or scoped credentials—if an agent can access them, so can an attacker.

### [​](#handling-secrets-safely)Handling secrets safely

If your agent needs to call authenticated APIs or access protected resources, you have two options:

{' ' * (self.list_depth - 1)}- 
**Keep secrets in tools outside the sandbox.** Define tools that run in your host environment (not inside the sandbox) and handle authentication there. The agent calls these tools by name, but never sees the credentials. This is the recommended approach.

{' ' * (self.list_depth - 1)}- 
**Use a network proxy that injects credentials.** Some sandbox providers support proxies that intercept outgoing HTTP requests from the sandbox and attach credentials (e.g., `Authorization` headers) before forwarding them. The agent never sees the secret—it just makes plain requests to a URL. This approach is not yet widely available across providers.

If you must inject secrets into a sandbox (not recommended), take these precautions:

{' ' * (self.list_depth - 1)}- Enable [human-in-the-loop](/oss/python/deepagents/human-in-the-loop) approval for **all** tool calls, not just sensitive ones

{' ' * (self.list_depth - 1)}- Block or restrict network access from the sandbox to limit exfiltration paths

{' ' * (self.list_depth - 1)}- Use the narrowest possible credential scope and shortest possible lifetime

{' ' * (self.list_depth - 1)}- Monitor sandbox network traffic for unexpected outbound requests
Even with these safeguards, this remains an unsafe workaround. A sufficiently creative enough context injection attack can bypass output filtering and HITL review.

### [​](#general-best-practices)General best practices

{' ' * (self.list_depth - 1)}- Review sandbox outputs before acting on them in your application

{' ' * (self.list_depth - 1)}- Block sandbox network access when not needed

{' ' * (self.list_depth - 1)}- Use [middleware](/oss/python/langchain/middleware) to filter or redact sensitive patterns in tool outputs

{' ' * (self.list_depth - 1)}- Treat everything produced inside the sandbox as untrusted input

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/sandboxes.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[SkillsPrevious](/oss/python/deepagents/skills)[InterpretersNext](/oss/python/deepagents/interpreters)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
