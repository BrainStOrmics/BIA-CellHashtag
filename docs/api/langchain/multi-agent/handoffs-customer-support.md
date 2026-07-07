# Build customer support with handoffs

[Tutorials](/oss/python/deepagents/data-analysis)[Multi-agent](/oss/python/langchain/multi-agent/subagents-personal-assistant)

# Build customer support with handoffs
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.The [state machine pattern](/oss/python/langchain/multi-agent/handoffs) describes workflows where an agent’s behavior changes as it moves through different states of a task. This tutorial shows how to implement a state machine by using tool calls to dynamically change a single agent’s configuration—updating its available tools and instructions based on the current state. The state can be determined from multiple sources: the agent’s past actions (tool calls), external state (such as API call results), or even initial user input (for example, by running a classifier to determine user intent).
In this tutorial, you’ll build a customer support agent that does the following:

{' ' * (self.list_depth - 1)}- Collects warranty information before proceeding.

{' ' * (self.list_depth - 1)}- Classifies issues as hardware or software.

{' ' * (self.list_depth - 1)}- Provides solutions or escalates to human support.

{' ' * (self.list_depth - 1)}- Maintains conversation state across multiple turns.

Unlike the [subagents pattern](/oss/python/langchain/multi-agent/subagents-personal-assistant) where sub-agents are called as tools, the **state machine pattern** uses a single agent whose configuration changes based on workflow progress. Each “step” is just a different configuration (system prompt + tools) of the same underlying agent, selected dynamically based on state.
Here’s the workflow we’ll build:

## [​](#setup)Setup

### [​](#installation)Installation

This tutorial requires the `langchain` package:
pipuvconda
```
pip install langchain

```

For more details, see our [Installation guide](/oss/python/langchain/install).

### [​](#langsmith)LangSmith

Set up [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-multi-agent-handoffs-customer-support) to inspect what is happening inside your agent. Then set the following environment variables:
bashpython
```
export LANGSMITH_TRACING="true"
export LANGSMITH_API_KEY="..."

```

### [​](#select-an-llm)Select an LLM

Select a chat model from LangChain’s suite of integrations:

{' ' * (self.list_depth - 1)}- OpenAI
{' ' * (self.list_depth - 1)}- Anthropic
{' ' * (self.list_depth - 1)}- Azure
{' ' * (self.list_depth - 1)}- Google Gemini
{' ' * (self.list_depth - 1)}- AWS Bedrock
{' ' * (self.list_depth - 1)}- HuggingFace
{' ' * (self.list_depth - 1)}- OpenRouter👉 Read the [OpenAI chat model integration docs](/oss/python/integrations/chat/openai)
```
pip install -U "langchain[openai]"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["OPENAI_API_KEY"] = "sk-..."

model = init_chat_model("gpt-5.4")

```
👉 Read the [Anthropic chat model integration docs](/oss/python/integrations/chat/anthropic)
```
pip install -U "langchain[anthropic]"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["ANTHROPIC_API_KEY"] = "sk-..."

model = init_chat_model("claude-sonnet-4-6")

```
👉 Read the [Azure chat model integration docs](/oss/python/integrations/chat/azure_chat_openai)
```
pip install -U "langchain[openai]"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["AZURE_OPENAI_API_KEY"] = "..."
os.environ["AZURE_OPENAI_ENDPOINT"] = "..."
os.environ["OPENAI_API_VERSION"] = "2025-03-01-preview"

model = init_chat_model(
 "azure_openai:gpt-5.4",
 azure_deployment=os.environ["AZURE_OPENAI_DEPLOYMENT_NAME"],
)

```
👉 Read the [Google GenAI chat model integration docs](/oss/python/integrations/chat/google_generative_ai)
```
pip install -U "langchain[google-genai]"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["GOOGLE_API_KEY"] = "..."

model = init_chat_model("google_genai:gemini-2.5-flash-lite")

```
👉 Read the [AWS Bedrock chat model integration docs](/oss/python/integrations/chat/bedrock)
```
pip install -U "langchain[aws]"

```
init_chat_modelModel Class
```
from langchain.chat_models import init_chat_model

# Follow the steps here to configure your credentials:
# https://docs.aws.amazon.com/bedrock/latest/userguide/getting-started.html

model = init_chat_model(
 "anthropic.claude-3-5-sonnet-20240620-v1:0",
 model_provider="bedrock_converse",
)

```
👉 Read the [HuggingFace chat model integration docs](/oss/python/integrations/chat/huggingface)
```
pip install -U "langchain[huggingface]"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["HUGGINGFACEHUB_API_TOKEN"] = "hf_..."

model = init_chat_model(
 "microsoft/Phi-3-mini-4k-instruct",
 model_provider="huggingface",
 temperature=0.7,
 max_tokens=1024,
)

```
👉 Read the [OpenRouter chat model integration docs](/oss/python/integrations/chat/openrouter)
```
pip install -U "langchain-openrouter"

```
init_chat_modelModel Class
```
import os
from langchain.chat_models import init_chat_model

os.environ["OPENROUTER_API_KEY"] = "sk-..."

model = init_chat_model(
 "auto",
 model_provider="openrouter",
)

```

## [​](#1-define-custom-state)1. Define custom state

First, define a custom state schema that tracks which step is currently active:

```
from langchain.agents import AgentState
from typing_extensions import NotRequired
from typing import Literal

# Define the possible workflow steps
SupportStep = Literal["warranty_collector", "issue_classifier", "resolution_specialist"]

class SupportState(AgentState):
 """State for customer support workflow."""
 current_step: NotRequired[SupportStep]
 warranty_status: NotRequired[Literal["in_warranty", "out_of_warranty"]]
 issue_type: NotRequired[Literal["hardware", "software"]]

```

The `current_step` field is the core of the state machine pattern - it determines which configuration (prompt + tools) is loaded on each turn.

## [​](#2-create-tools-that-manage-workflow-state)2. Create tools that manage workflow state

Create tools that update the workflow state. These tools allow the agent to record information and transition to the next step.
The key is using `Command` to update state, including the `current_step` field:

```
from langchain.tools import tool, ToolRuntime
from langchain.messages import ToolMessage
from langgraph.types import Command

@tool
def record_warranty_status(
 status: Literal["in_warranty", "out_of_warranty"],
 runtime: ToolRuntime[None, SupportState],
) -> Command:
 """Record the customer's warranty status and transition to issue classification."""
 return Command(
 update={
 "messages": [
 ToolMessage(
 content=f"Warranty status recorded as: {status}",
 tool_call_id=runtime.tool_call_id,
 )
 ],
 "warranty_status": status,
 "current_step": "issue_classifier",
 }
 )

@tool
def record_issue_type(
 issue_type: Literal["hardware", "software"],
 runtime: ToolRuntime[None, SupportState],
) -> Command:
 """Record the type of issue and transition to resolution specialist."""
 return Command(
 update={
 "messages": [
 ToolMessage(
 content=f"Issue type recorded as: {issue_type}",
 tool_call_id=runtime.tool_call_id,
 )
 ],
 "issue_type": issue_type,
 "current_step": "resolution_specialist",
 }
 )

@tool
def escalate_to_human(reason: str) -> str:
 """Escalate the case to a human support specialist."""
 # In a real system, this would create a ticket, notify staff, etc.
 return f"Escalating to human support. Reason: {reason}"

@tool
def provide_solution(solution: str) -> str:
 """Provide a solution to the customer's issue."""
 return f"Solution provided: {solution}"

```

Notice how `record_warranty_status` and `record_issue_type` return `Command` objects that update both the data (`warranty_status`, `issue_type`) AND the `current_step`. This is how the state machine works - tools control workflow progression.

## [​](#3-define-step-configurations)3. Define step configurations

Define prompts and tools for each step. First, define the prompts for each step:

View complete prompt definitions
```
# Define prompts as constants for easy reference
WARRANTY_COLLECTOR_PROMPT = """You are a customer support agent helping with device issues.

CURRENT STAGE: Warranty verification

At this step, you need to:
1. Greet the customer warmly
2. Ask if their device is under warranty
3. Use record_warranty_status to record their response and move to the next step

Be conversational and friendly. Don't ask multiple questions at once."""

ISSUE_CLASSIFIER_PROMPT = """You are a customer support agent helping with device issues.

CURRENT STAGE: Issue classification
CUSTOMER INFO: Warranty status is {warranty_status}

At this step, you need to:
1. Ask the customer to describe their issue
2. Determine if it's a hardware issue (physical damage, broken parts) or software issue (app crashes, performance)
3. Use record_issue_type to record the classification and move to the next step

If unclear, ask clarifying questions before classifying."""

RESOLUTION_SPECIALIST_PROMPT = """You are a customer support agent helping with device issues.

CURRENT STAGE: Resolution
CUSTOMER INFO: Warranty status is {warranty_status}, issue type is {issue_type}

At this step, you need to:
1. For SOFTWARE issues: provide troubleshooting steps using provide_solution
2. For HARDWARE issues:
 - If IN WARRANTY: explain warranty repair process using provide_solution
 - If OUT OF WARRANTY: escalate_to_human for paid repair options

Be specific and helpful in your solutions."""

```

Then map step names to their configurations using a dictionary:

```
# Step configuration: maps step name to (prompt, tools, required_state)
STEP_CONFIG = {
 "warranty_collector": {
 "prompt": WARRANTY_COLLECTOR_PROMPT,
 "tools": [record_warranty_status],
 "requires": [],
 },
 "issue_classifier": {
 "prompt": ISSUE_CLASSIFIER_PROMPT,
 "tools": [record_issue_type],
 "requires": ["warranty_status"],
 },
 "resolution_specialist": {
 "prompt": RESOLUTION_SPECIALIST_PROMPT,
 "tools": [provide_solution, escalate_to_human],
 "requires": ["warranty_status", "issue_type"],
 },
}

```

This dictionary-based configuration makes it easy to:

{' ' * (self.list_depth - 1)}- See all steps at a glance

{' ' * (self.list_depth - 1)}- Add new steps (just add another entry)

{' ' * (self.list_depth - 1)}- Understand the workflow dependencies (`requires` field)

{' ' * (self.list_depth - 1)}- Use prompt templates with state variables (e.g., `{warranty_status}`)

## [​](#4-create-step-based-middleware)4. Create step-based middleware

Create middleware that reads `current_step` from state and applies the appropriate configuration. We’ll use the `@wrap_model_call` decorator for a clean implementation:

```
from langchain.agents.middleware import wrap_model_call, ModelRequest, ModelResponse
from typing import Callable

@wrap_model_call
def apply_step_config(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
) -> ModelResponse:
 """Configure agent behavior based on the current step."""
 # Get current step (defaults to warranty_collector for first interaction)
 current_step = request.state.get("current_step", "warranty_collector")

 # Look up step configuration
 stage_config = STEP_CONFIG[current_step]

 # Validate required state exists
 for key in stage_config["requires"]:
 if request.state.get(key) is None:
 raise ValueError(f"{key} must be set before reaching {current_step}")

 # Format prompt with state values (supports {warranty_status}, {issue_type}, etc.)
 system_prompt = stage_config["prompt"].format(**request.state)

 # Inject system prompt and step-specific tools
 request = request.override(
 system_prompt=system_prompt,
 tools=stage_config["tools"],
 )

 return handler(request)

```

This middleware:

{' ' * (self.list_depth - 1)}- **Reads current step**: Gets `current_step` from state (defaults to “warranty_collector”).

{' ' * (self.list_depth - 1)}- **Looks up configuration**: Finds the matching entry in `STEP_CONFIG`.

{' ' * (self.list_depth - 1)}- **Validates dependencies**: Ensures required state fields exist.

{' ' * (self.list_depth - 1)}- **Formats prompt**: Injects state values into the prompt template.

{' ' * (self.list_depth - 1)}- **Applies configuration**: Overrides the system prompt and available tools.

The `request.override()` method is key - it allows us to dynamically change the agent’s behavior based on state without creating separate agent instances.

## [​](#5-create-the-agent)5. Create the agent

Now create the agent with the step-based middleware and a checkpointer for state persistence:

```
from langchain.agents import create_agent
from langgraph.checkpoint.memory import InMemorySaver

# Collect all tools from all step configurations
all_tools = [
 record_warranty_status,
 record_issue_type,
 provide_solution,
 escalate_to_human,
]

# Create the agent with step-based configuration
agent = create_agent(
 model,
 tools=all_tools,
 state_schema=SupportState,
 middleware=[apply_step_config],
 checkpointer=InMemorySaver(),
)

```

**Why a checkpointer?** The checkpointer maintains state across conversation turns. Without it, the `current_step` state would be lost between user messages, breaking the workflow.

## [​](#6-test-the-workflow)6. Test the workflow

Test the complete workflow:

```
from langchain.messages import HumanMessage
from langchain_core.utils.uuid import uuid7

# Configuration for this conversation thread
thread_id = str(uuid7())
config = {"configurable": {"thread_id": thread_id}}

# Turn 1: Initial message - starts with warranty_collector step
print("=== Turn 1: Warranty Collection ===")
result = agent.invoke(
 {"messages": [HumanMessage("Hi, my phone screen is cracked")]},
 config
)
for msg in result['messages']:
 msg.pretty_print()

# Turn 2: User responds about warranty
print("\n=== Turn 2: Warranty Response ===")
result = agent.invoke(
 {"messages": [HumanMessage("Yes, it's still under warranty")]},
 config
)
for msg in result['messages']:
 msg.pretty_print()
print(f"Current step: {result.get('current_step')}")

# Turn 3: User describes the issue
print("\n=== Turn 3: Issue Description ===")
result = agent.invoke(
 {"messages": [HumanMessage("The screen is physically cracked from dropping it")]},
 config
)
for msg in result['messages']:
 msg.pretty_print()
print(f"Current step: {result.get('current_step')}")

# Turn 4: Resolution
print("\n=== Turn 4: Resolution ===")
result = agent.invoke(
 {"messages": [HumanMessage("What should I do?")]},
 config
)
for msg in result['messages']:
 msg.pretty_print()

```

Expected flow:

{' ' * (self.list_depth - 1)}- **Warranty verification step**: Asks about warranty status

{' ' * (self.list_depth - 1)}- **Issue classification step**: Asks about the problem, determines it’s hardware

{' ' * (self.list_depth - 1)}- **Resolution step**: Provides warranty repair instructions

## [​](#7-understanding-state-transitions)7. Understanding state transitions

Let’s trace what happens at each turn:

### [​](#turn-1-initial-message)Turn 1: Initial message

```
{
 "messages": [HumanMessage("Hi, my phone screen is cracked")],
 "current_step": "warranty_collector" # Default value
}

```

Middleware applies:

{' ' * (self.list_depth - 1)}- System prompt: `WARRANTY_COLLECTOR_PROMPT`

{' ' * (self.list_depth - 1)}- Tools: `[record_warranty_status]`

### [​](#turn-2-after-warranty-recorded)Turn 2: After warranty recorded

Tool call: `record_warranty_status("in_warranty")` returns:

```
Command(update={
 "warranty_status": "in_warranty",
 "current_step": "issue_classifier" # State transition!
})

```

Next turn, middleware applies:

{' ' * (self.list_depth - 1)}- System prompt: `ISSUE_CLASSIFIER_PROMPT` (formatted with `warranty_status="in_warranty"`)

{' ' * (self.list_depth - 1)}- Tools: `[record_issue_type]`

### [​](#turn-3-after-issue-classified)Turn 3: After issue classified

Tool call: `record_issue_type("hardware")` returns:

```
Command(update={
 "issue_type": "hardware",
 "current_step": "resolution_specialist" # State transition!
})

```

Next turn, middleware applies:

{' ' * (self.list_depth - 1)}- System prompt: `RESOLUTION_SPECIALIST_PROMPT` (formatted with `warranty_status` and `issue_type`)

{' ' * (self.list_depth - 1)}- Tools: `[provide_solution, escalate_to_human]`

The key insight: **Tools drive the workflow** by updating `current_step`, and **middleware responds** by applying the appropriate configuration on the next turn.

## [​](#8-manage-message-history)8. Manage message history

As the agent progresses through steps, message history grows. Use [summarization middleware](/oss/python/langchain/short-term-memory#summarize-messages) to compress earlier messages while preserving conversational context:

```
from langchain.agents import create_agent
from langchain.agents.middleware import SummarizationMiddleware 
from langgraph.checkpoint.memory import InMemorySaver

agent = create_agent(
 model,
 tools=all_tools,
 state_schema=SupportState,
 middleware=[
 apply_step_config,
 SummarizationMiddleware(
 model="gpt-5.4-mini",
 trigger=("tokens", 4000),
 keep=("messages", 10)
 )
 ],
 checkpointer=InMemorySaver(),
)

```

See the [short-term memory guide](/oss/python/langchain/short-term-memory) for other memory management techniques.

## [​](#9-add-flexibility-go-back)9. Add flexibility: Go back

Some workflows need to allow users to return to previous steps to correct information (e.g., changing warranty status or issue classification). However, not all transitions make sense—for example, you typically can’t go back once a refund has been processed. For this support workflow, we’ll add tools to return to the warranty verification and issue classification steps.
If your workflow requires arbitrary transitions between most steps, consider whether you need a structured workflow at all. This pattern works best when steps follow a clear sequential progression with occasional backwards transitions for corrections.
Add “go back” tools to the resolution step:

```
@tool
def go_back_to_warranty() -> Command:
 """Go back to warranty verification step."""
 return Command(update={"current_step": "warranty_collector"})

@tool
def go_back_to_classification() -> Command:
 """Go back to issue classification step."""
 return Command(update={"current_step": "issue_classifier"})

# Update the resolution_specialist configuration to include these tools
STEP_CONFIG["resolution_specialist"]["tools"].extend([
 go_back_to_warranty,
 go_back_to_classification
])

```

Update the resolution specialist’s prompt to mention these tools:

```
RESOLUTION_SPECIALIST_PROMPT = """You are a customer support agent helping with device issues.

CURRENT STAGE: Resolution
CUSTOMER INFO: Warranty status is {warranty_status}, issue type is {issue_type}

At this step, you need to:
1. For SOFTWARE issues: provide troubleshooting steps using provide_solution
2. For HARDWARE issues:
 - If IN WARRANTY: explain warranty repair process using provide_solution
 - If OUT OF WARRANTY: escalate_to_human for paid repair options

If the customer indicates any information was wrong, use:
- go_back_to_warranty to correct warranty status
- go_back_to_classification to correct issue type

Be specific and helpful in your solutions."""

```

Now the agent can handle corrections:

```
result = agent.invoke(
 {"messages": [HumanMessage("Actually, I made a mistake - my device is out of warranty")]},
 config
)
# Agent will call go_back_to_warranty and restart the warranty verification step

```

## [​](#complete-example)Complete example

Here’s everything together in a runnable script:

Show Complete code
```
"""
Customer Support State Machine Example

This example demonstrates the state machine pattern.
A single agent dynamically changes its behavior based on the current_step state,
creating a state machine for sequential information collection.
"""

from langchain_core.utils.uuid import uuid7

from langgraph.checkpoint.memory import InMemorySaver
from langgraph.types import Command
from typing import Callable, Literal
from typing_extensions import NotRequired

from langchain.agents import AgentState, create_agent
from langchain.agents.middleware import wrap_model_call, ModelRequest, ModelResponse, SummarizationMiddleware
from langchain.chat_models import init_chat_model
from langchain.messages import HumanMessage, ToolMessage
from langchain.tools import tool, ToolRuntime

model = init_chat_model("google_genai:gemini-3.5-flash")

# Define the possible workflow steps
SupportStep = Literal["warranty_collector", "issue_classifier", "resolution_specialist"]

class SupportState(AgentState):
 """State for customer support workflow."""

 current_step: NotRequired[SupportStep]
 warranty_status: NotRequired[Literal["in_warranty", "out_of_warranty"]]
 issue_type: NotRequired[Literal["hardware", "software"]]

@tool
def record_warranty_status(
 status: Literal["in_warranty", "out_of_warranty"],
 runtime: ToolRuntime[None, SupportState],
) -> Command:
 """Record the customer's warranty status and transition to issue classification."""
 return Command(
 update={
 "messages": [
 ToolMessage(
 content=f"Warranty status recorded as: {status}",
 tool_call_id=runtime.tool_call_id,
 )
 ],
 "warranty_status": status,
 "current_step": "issue_classifier",
 }
 )

@tool
def record_issue_type(
 issue_type: Literal["hardware", "software"],
 runtime: ToolRuntime[None, SupportState],
) -> Command:
 """Record the type of issue and transition to resolution specialist."""
 return Command(
 update={
 "messages": [
 ToolMessage(
 content=f"Issue type recorded as: {issue_type}",
 tool_call_id=runtime.tool_call_id,
 )
 ],
 "issue_type": issue_type,
 "current_step": "resolution_specialist",
 }
 )

@tool
def escalate_to_human(reason: str) -> str:
 """Escalate the case to a human support specialist."""
 # In a real system, this would create a ticket, notify staff, etc.
 return f"Escalating to human support. Reason: {reason}"

@tool
def provide_solution(solution: str) -> str:
 """Provide a solution to the customer's issue."""
 return f"Solution provided: {solution}"

# Define prompts as constants
WARRANTY_COLLECTOR_PROMPT = """You are a customer support agent helping with device issues.

CURRENT STEP: Warranty verification

At this step, you need to:
1. Greet the customer warmly
2. Ask if their device is under warranty
3. Use record_warranty_status to record their response and move to the next step

Be conversational and friendly. Don't ask multiple questions at once."""

ISSUE_CLASSIFIER_PROMPT = """You are a customer support agent helping with device issues.

CURRENT STEP: Issue classification
CUSTOMER INFO: Warranty status is {warranty_status}

At this step, you need to:
1. Ask the customer to describe their issue
2. Determine if it's a hardware issue (physical damage, broken parts) or software issue (app crashes, performance)
3. Use record_issue_type to record the classification and move to the next step

If unclear, ask clarifying questions before classifying."""

RESOLUTION_SPECIALIST_PROMPT = """You are a customer support agent helping with device issues.

CURRENT STEP: Resolution
CUSTOMER INFO: Warranty status is {warranty_status}, issue type is {issue_type}

At this step, you need to:
1. For SOFTWARE issues: provide troubleshooting steps using provide_solution
2. For HARDWARE issues:
 - If IN WARRANTY: explain warranty repair process using provide_solution
 - If OUT OF WARRANTY: escalate_to_human for paid repair options

Be specific and helpful in your solutions."""

# Step configuration: maps step name to (prompt, tools, required_state)
STEP_CONFIG = {
 "warranty_collector": {
 "prompt": WARRANTY_COLLECTOR_PROMPT,
 "tools": [record_warranty_status],
 "requires": [],
 },
 "issue_classifier": {
 "prompt": ISSUE_CLASSIFIER_PROMPT,
 "tools": [record_issue_type],
 "requires": ["warranty_status"],
 },
 "resolution_specialist": {
 "prompt": RESOLUTION_SPECIALIST_PROMPT,
 "tools": [provide_solution, escalate_to_human],
 "requires": ["warranty_status", "issue_type"],
 },
}

@wrap_model_call
def apply_step_config(
 request: ModelRequest,
 handler: Callable[[ModelRequest], ModelResponse],
) -> ModelResponse:
 """Configure agent behavior based on the current step."""
 # Get current step (defaults to warranty_collector for first interaction)
 current_step = request.state.get("current_step", "warranty_collector")

 # Look up step configuration
 step_config = STEP_CONFIG[current_step]

 # Validate required state exists
 for key in step_config["requires"]:
 if request.state.get(key) is None:
 raise ValueError(f"{key} must be set before reaching {current_step}")

 # Format prompt with state values
 system_prompt = step_config["prompt"].format(**request.state)

 # Inject system prompt and step-specific tools
 request = request.override(
 system_prompt=system_prompt,
 tools=step_config["tools"],
 )

 return handler(request)

# Collect all tools from all step configurations
all_tools = [
 record_warranty_status,
 record_issue_type,
 provide_solution,
 escalate_to_human,
]

# Create the agent with step-based configuration and summarization
agent = create_agent(
 model,
 tools=all_tools,
 state_schema=SupportState,
 middleware=[
 apply_step_config,
 SummarizationMiddleware(
 model="gpt-5.4-mini",
 trigger=("tokens", 4000),
 keep=("messages", 10)
 )
 ],
 checkpointer=InMemorySaver(),
)

# ============================================================================
# Test the workflow
# ============================================================================

if __name__ == "__main__":
 thread_id = str(uuid7())
 config = {"configurable": {"thread_id": thread_id}}

 result = agent.invoke(
 {"messages": [HumanMessage("Hi, my phone screen is cracked")]},
 config
 )

 result = agent.invoke(
 {"messages": [HumanMessage("Yes, it's still under warranty")]},
 config
 )

 result = agent.invoke(
 {"messages": [HumanMessage("The screen is physically cracked from dropping it")]},
 config
 )

 result = agent.invoke(
 {"messages": [HumanMessage("What should I do?")]},
 config
 )
 for msg in result['messages']:
 msg.pretty_print()

```

## [​](#next-steps)Next steps

{' ' * (self.list_depth - 1)}- Learn about the [subagents pattern](/oss/python/langchain/multi-agent/subagents-personal-assistant) for centralized orchestration

{' ' * (self.list_depth - 1)}- Explore [middleware](/oss/python/langchain/middleware) for more dynamic behaviors

{' ' * (self.list_depth - 1)}- Read the [multi-agent overview](/oss/python/langchain/multi-agent) to compare patterns

{' ' * (self.list_depth - 1)}- Use [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-multi-agent-handoffs-customer-support) to debug and monitor your multi-agent system

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/multi-agent/handoffs-customer-support.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Build a personal assistant with subagentsPrevious](/oss/python/langchain/multi-agent/subagents-personal-assistant)[Build a multi-source knowledge base with routingNext](/oss/python/langchain/multi-agent/router-knowledge-base)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
