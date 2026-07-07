# Human-in-the-loop

[Core capabilities](/oss/python/deepagents/harness)

# Human-in-the-loop
Copy page

Learn how to configure human approval for sensitive tool operationsCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Some tool operations may be sensitive and require human approval before execution. Deep Agents support human-in-the-loop workflows through LangGraph’s interrupt capabilities. You can configure which tools require approval using the `interrupt_on` parameter.

## [​](#basic-configuration)Basic configuration

The `interrupt_on` parameter accepts a dictionary mapping tool names to interrupt configurations. Each tool can be configured with:

{' ' * (self.list_depth - 1)}- **`True`**: Enable interrupts with default behavior (approve, edit, reject, respond allowed)

{' ' * (self.list_depth - 1)}- **`False`**: Disable interrupts for this tool

{' ' * (self.list_depth - 1)}- **`{"allowed_decisions": [...]}`**: Custom configuration with specific allowed decisions

```
from langchain.tools import tool
from deepagents import create_deep_agent
from langgraph.checkpoint.memory import MemorySaver

@tool
def remove_file(path: str) -> str:
 """Delete a file from the filesystem."""
 return f"Deleted {path}"

@tool
def fetch_file(path: str) -> str:
 """Read a file from the filesystem."""
 return f"Contents of {path}"

@tool
def notify_email(to: str, subject: str, body: str) -> str:
 """Send an email."""
 return f"Sent email to {to}"

# Checkpointer is REQUIRED for human-in-the-loop
checkpointer = MemorySaver()

agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[remove_file, fetch_file, notify_email],
 interrupt_on={
 "remove_file": True, # Default: approve, edit, reject, respond
 "fetch_file": False, # No interrupts needed
 "notify_email": {"allowed_decisions": ["approve", "reject"]}, # No editing
 },
 checkpointer=checkpointer, # Required!
)

```

## [​](#decision-types)Decision types

The `allowed_decisions` list controls what actions a human can take when reviewing a tool call:

{' ' * (self.list_depth - 1)}- **`"approve"`**: Execute the tool with the original arguments as proposed by the agent

{' ' * (self.list_depth - 1)}- **`"edit"`**: Modify the tool arguments before execution

{' ' * (self.list_depth - 1)}- **`"reject"`**: Skip executing this tool call entirely

{' ' * (self.list_depth - 1)}- **`"respond"`**: Return the human’s message directly as the tool result, skipping execution — for “ask user” style tools

You can customize which decisions are available for each tool:

```
interrupt_on = {
 # Sensitive operations: allow all options
 "delete_file": {"allowed_decisions": ["approve", "edit", "reject"]},

 # Moderate risk: approval or rejection only
 "write_file": {"allowed_decisions": ["approve", "reject"]},

 # Must approve (no rejection allowed)
 "critical_operation": {"allowed_decisions": ["approve"]},
}

```

## [​](#handle-interrupts)Handle interrupts

When an interrupt is triggered, the agent pauses execution and returns control. Check for interrupts in the result and handle them accordingly.

```
from langchain_core.utils.uuid import uuid7
from langgraph.types import Command

# Create config with thread_id for state persistence
config = {"configurable": {"thread_id": str(uuid7())}}

# Invoke the agent
result = agent.invoke(
 {"messages": [{"role": "user", "content": "Delete the file temp.txt"}]},
 config=config,
 version="v2",
)

# Check if execution was interrupted
if result.interrupts:
 # Extract interrupt information
 interrupt_value = result.interrupts[0].value 
 action_requests = interrupt_value["action_requests"]
 review_configs = interrupt_value["review_configs"]

 # Create a lookup map from tool name to review config
 config_map = {cfg["action_name"]: cfg for cfg in review_configs}

 # Display the pending actions to the user
 for action in action_requests:
 review_config = config_map[action["name"]]
 print(f"Tool: {action['name']}")
 print(f"Arguments: {action['args']}")
 print(f"Allowed decisions: {review_config['allowed_decisions']}")

 # Get user decisions (one per action_request, in order)
 decisions = [
 {"type": "approve"} # User approved the deletion
 ]

 # Resume execution with decisions
 result = agent.invoke(
 Command(resume={"decisions": decisions}),
 config=config, # Must use the same config!
 version="v2",
 )

# Process final result
print(result.value["messages"][-1].content)

```

## [​](#multiple-tool-calls)Multiple tool calls

When the agent calls multiple tools that require approval, all interrupts are batched together in a single interrupt. You must provide decisions for each one in order.

```
config = {"configurable": {"thread_id": str(uuid7())}}

result = agent.invoke(
 {"messages": [{
 "role": "user",
 "content": "Delete temp.txt and send an email to admin@example.com"
 }]},
 config=config,
 version="v2",
)

if result.interrupts:
 interrupt_value = result.interrupts[0].value 
 action_requests = interrupt_value["action_requests"]

 # Two tools need approval
 assert len(action_requests) == 2

 # Provide decisions in the same order as action_requests
 decisions = [
 {"type": "approve"}, # First tool: delete_file
 {"type": "reject"} # Second tool: send_email
 ]

 result = agent.invoke(
 Command(resume={"decisions": decisions}),
 config=config,
 version="v2",
 )

```

## [​](#edit-tool-arguments)Edit tool arguments

When `"edit"` is in the allowed decisions, you can modify the tool arguments before execution:

```
if result.interrupts:
 interrupt_value = result.interrupts[0].value 
 action_request = interrupt_value["action_requests"][0]

 # Original args from the agent
 print(action_request["args"]) # {"to": "everyone@company.com", ...}

 # User decides to edit the recipient
 decisions = [{
 "type": "edit",
 "edited_action": {
 "name": action_request["name"], # Must include the tool name
 "args": {"to": "team@company.com", "subject": "...", "body": "..."}
 }
 }]

 result = agent.invoke(
 Command(resume={"decisions": decisions}),
 config=config,
 version="v2",
 )

```

## [​](#subagent-interrupts)Subagent interrupts

When using subagents, you can use interrupts [on tool calls](#interrupts-on-tool-calls) and [within tool calls](#interrupts-within-tool-calls).

### [​](#interrupts-on-tool-calls)Interrupts on tool calls

Each subagent can have its own `interrupt_on` configuration that overrides the main agent’s settings:

```
agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[delete_file, read_file],
 interrupt_on={
 "delete_file": True,
 "read_file": False,
 },
 subagents=[{
 "name": "file-manager",
 "description": "Manages file operations",
 "system_prompt": "You are a file management assistant.",
 "tools": [delete_file, read_file],
 "interrupt_on": {
 # Override: require approval for reads in this subagent
 "delete_file": True,
 "read_file": True, # Different from main agent!
 }
 }],
 checkpointer=checkpointer
)

```

When a subagent triggers an interrupt, the handling is the same—check for `interrupts` on the result and resume with `Command`.

### [​](#interrupts-within-tool-calls)Interrupts within tool calls

Subagent tools can call `interrupt()` directly to pause execution and await approval:

```
from langchain.agents import create_agent
from langchain_anthropic import ChatAnthropic
from langchain.messages import HumanMessage
from langchain.tools import tool
from langgraph.checkpoint.memory import InMemorySaver
from langgraph.types import Command, interrupt

from deepagents.graph import create_deep_agent
from deepagents.middleware.subagents import CompiledSubAgent

@tool(description="Request human approval before proceeding with an action.")
def request_approval(action_description: str) -> str:
 """Request human approval using the interrupt() primitive."""
 # interrupt() pauses execution and returns the value passed to Command(resume=...)
 approval = interrupt({
 "type": "approval_request",
 "action": action_description,
 "message": f"Please approve or reject: {action_description}",
 })

 if approval.get("approved"):
 return f"Action '{action_description}' was APPROVED. Proceeding..."
 else:
 return f"Action '{action_description}' was REJECTED. Reason: {approval.get('reason', 'No reason provided')}"

def main():
 checkpointer = InMemorySaver()
 model = ChatAnthropic(
 model_name="claude-sonnet-4-6",
 max_tokens=4096,
 )

 compiled_subagent = create_agent(
 model=model,
 tools=[request_approval],
 name="approval-agent",
 )

 parent_agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 checkpointer=checkpointer,
 subagents=[
 CompiledSubAgent(
 name="approval-agent",
 description="An agent that can request approvals",
 runnable=compiled_subagent,
 )
 ],
 )

 thread_id = "test_interrupt_directly"
 config = {"configurable": {"thread_id": thread_id}}

 print("Invoking agent - sub-agent will use request_approval tool...")

 result = parent_agent.invoke(
 {
 "messages": [
 HumanMessage(
 content="Use the task tool to launch the approval-agent sub-agent. "
 "Tell it to use the request_approval tool to request approval for 'deploying to production'."
 )
 ]
 },
 config=config,
 version="v2",
 )

 # Check for interrupt
 if result.interrupts:
 interrupt_value = result.interrupts[0].value 
 print(f"\nInterrupt received!")
 print(f" Type: {interrupt_value.get('type')}")
 print(f" Action: {interrupt_value.get('action')}")
 print(f" Message: {interrupt_value.get('message')}")

 print("\nResuming with Command(resume={'approved': True})...")
 result2 = parent_agent.invoke(
 Command(resume={"approved": True}),
 config=config,
 version="v2",
 )

 if not result2.interrupts:
 print("\nExecution completed!")
 # Find the tool response
 tool_msgs = [m for m in result2.value.get("messages", []) if m.type == "tool"]
 if tool_msgs:
 print(f" Tool result: {tool_msgs[-1].content}")
 else:
 print("\nAnother interrupt occurred")
 else:
 print("\n No interrupt - the model may not have called request_approval")

if __name__ == "__main__":
 main()

```

When run, this produces the following output:

```
Invoking agent - sub-agent will use request_approval tool...

Interrupt received!
 Type: approval_request
 Action: deploying to production
 Message: Please approve or reject: deploying to production

Resuming with Command(resume={'approved': True})...

Execution completed!
 Tool result: Great! The approval request has been processed. The action **"deploying to production"** was **APPROVED**. You can now proceed with the production deployment.

```

## [​](#best-practices)Best practices

### [​](#always-use-a-checkpointer)Always use a checkpointer

Human-in-the-loop requires a checkpointer to persist agent state between the interrupt and resume:

```
from langgraph.checkpoint.memory import MemorySaver

checkpointer = MemorySaver()
agent = create_deep_agent(
 model="google_genai:gemini-3.5-flash",
 tools=[...],
 interrupt_on={...},
 checkpointer=checkpointer # Required for HITL
)

```

### [​](#use-the-same-thread-id)Use the same thread ID

When resuming, you must use the same config with the same `thread_id`:

```
# First call
config = {"configurable": {"thread_id": "my-thread"}}
result = agent.invoke(input, config=config, version="v2")

# Resume (use same config)
result = agent.invoke(Command(resume={...}), config=config, version="v2")

```

### [​](#match-decision-order-to-actions)Match decision order to actions

The decisions list must match the order of `action_requests`:

```
if result.interrupts:
 interrupt_value = result.interrupts[0].value 
 action_requests = interrupt_value["action_requests"]

 # Create one decision per action, in order
 decisions = []
 for action in action_requests:
 decision = get_user_decision(action) # Your logic
 decisions.append(decision)

 result = agent.invoke(
 Command(resume={"decisions": decisions}),
 config=config,
 version="v2",
 )

```

### [​](#tailor-configurations-by-risk)Tailor configurations by risk

Configure different tools based on their risk level:

```
interrupt_on = {
 # High risk: full control (approve, edit, reject)
 "delete_file": {"allowed_decisions": ["approve", "edit", "reject"]},
 "send_email": {"allowed_decisions": ["approve", "edit", "reject"]},

 # Medium risk: no editing allowed
 "write_file": {"allowed_decisions": ["approve", "reject"]},

 # Low risk: no interrupts
 "read_file": False,
 "list_files": False,
}

```

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/deepagents/human-in-the-loop.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Async subagentsPrevious](/oss/python/deepagents/async-subagents)[PermissionsNext](/oss/python/deepagents/permissions)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
