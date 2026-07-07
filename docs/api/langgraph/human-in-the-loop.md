# Interrupts

[Capabilities](/oss/python/langgraph/persistence)

# Interrupts
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.Interrupts allow you to pause graph execution at specific points and wait for external input before continuing. This enables human-in-the-loop patterns where you need external input to proceed. When an interrupt is triggered, LangGraph saves the graph state using its [persistence](/oss/python/langgraph/persistence) layer and waits indefinitely until you resume execution.
Interrupts work by calling the `interrupt()` function at any point in your graph nodes. The function accepts any JSON-serializable value which is surfaced to the caller. When you’re ready to continue, you resume execution by re-invoking the graph using `Command`, which then becomes the return value of the `interrupt()` call from inside the node.
Unlike static breakpoints (which pause before or after specific nodes), interrupts are **dynamic**: they can be placed anywhere in your code and can be conditional based on your application logic.

{' ' * (self.list_depth - 1)}- **Checkpointing keeps your place:** the checkpointer writes the exact graph state so you can resume later, even when in an error state.

{' ' * (self.list_depth - 1)}- **`thread_id` is your pointer:** set `config={"configurable": {"thread_id": ...}}` to tell the checkpointer which state to load.

{' ' * (self.list_depth - 1)}- **Interrupt payloads surface via `chunk["interrupts"]`:** when streaming with `version="v2"`, the values you pass to `interrupt()` appear in the `interrupts` field of `values` stream parts so you know what the graph is waiting on.

The `thread_id` you choose is effectively your persistent cursor. Reusing it resumes the same checkpoint; using a new value starts a brand-new thread with an empty state.

## [​](#pause-using-interrupt)Pause using `interrupt`

The [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) function pauses graph execution and returns a value to the caller. When you call [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) within a node, LangGraph saves the current graph state and waits for you to resume execution with input.
To use [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt), you need:

{' ' * (self.list_depth - 1)}- A **checkpointer** to persist the graph state (use a durable checkpointer in production)

{' ' * (self.list_depth - 1)}- A **thread ID** in your config so the runtime knows which state to resume from

{' ' * (self.list_depth - 1)}- To call `interrupt()` where you want to pause (payload must be JSON-serializable)

```
from langgraph.types import interrupt

def approval_node(state: State):
 # Pause and ask for approval
 approved = interrupt("Do you approve this action?")

 # When you resume, Command(resume=...) returns that value here
 return {"approved": approved}

```

When you call [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt), here’s what happens:

{' ' * (self.list_depth - 1)}- 
**Graph execution gets suspended** at the exact point where [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) is called

{' ' * (self.list_depth - 1)}- 
**State is saved** using the checkpointer so execution can be resumed later, In production, this should be a persistent checkpointer (e.g. backed by a database)

{' ' * (self.list_depth - 1)}- 
**Value is returned** to the caller on `result.interrupts` when using `version="v2"` (or under `__interrupt__` with the default v1 invoke API); it can be any JSON-serializable value (string, object, array, etc.)

{' ' * (self.list_depth - 1)}- 
**Graph waits indefinitely** until you resume execution with a response

{' ' * (self.list_depth - 1)}- 
**Response is passed back** into the node when you resume, becoming the return value of the `interrupt()` call

## [​](#resuming-interrupts)Resuming interrupts

After an interrupt pauses execution, you resume the graph by invoking it again with a `Command` that contains the resume value. The resume value is passed back to the `interrupt` call, allowing the node to continue execution with the external input.

{' ' * (self.list_depth - 1)}- v2 (LangGraph >= 1.1)
{' ' * (self.list_depth - 1)}- v1 (default)
```
from langgraph.types import Command

# Initial run - hits the interrupt and pauses
# thread_id is the persistent pointer (stores a stable ID in production)
config = {"configurable": {"thread_id": "thread-1"}}
result = graph.invoke({"input": "data"}, config=config, version="v2")

# result is a GraphOutput with .value and .interrupts
# .interrupts contains the payloads passed to interrupt()
print(result.interrupts)
# > (Interrupt(value='Do you approve this action?'),)

# Resume with the human's response
# The resume payload becomes the return value of interrupt() inside the node
graph.invoke(Command(resume=True), config=config, version="v2")

```

```
from langgraph.types import Command

config = {"configurable": {"thread_id": "thread-1"}}
result = graph.invoke({"input": "data"}, config=config)

# __interrupt__ contains the payload that was passed to interrupt()
print(result["__interrupt__"])
# > [Interrupt(value='Do you approve this action?')]

# Resume with the human's response
graph.invoke(Command(resume=True), config=config)

```

**Key points about resuming:**

{' ' * (self.list_depth - 1)}- You must use the **same thread ID** when resuming that was used when the interrupt occurred

{' ' * (self.list_depth - 1)}- The value passed to `Command(resume=...)` becomes the return value of the [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) call

{' ' * (self.list_depth - 1)}- The node restarts from the beginning of the node where the [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) was called when resumed, so any code before the [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) runs again

{' ' * (self.list_depth - 1)}- You can pass any JSON-serializable value as the resume value

`Command(resume=...)` is the **only** `Command` pattern intended as input to `invoke()`/`stream()`. The other `Command` parameters (`update`, `goto`, `graph`) are designed for [returning from node functions](/oss/python/langgraph/graph-api#command). Do not pass `Command(update=...)` as input to continue multi-turn conversations—pass a plain input dict instead.

## [​](#common-patterns)Common patterns

The key thing that interrupts unlock is the ability to pause execution and wait for external input. This is useful for a variety of use cases, including:

{' ' * (self.list_depth - 1)}- [Approval workflows](#approve-or-reject): Pause before executing critical actions (API calls, database changes, financial transactions)

{' ' * (self.list_depth - 1)}- [Handling multiple interrupts](#handling-multiple-interrupts): Pair interrupt IDs with resume values when resuming multiple interrupts in a single invocation

{' ' * (self.list_depth - 1)}- [Review and edit](#review-and-edit-state): Let humans review and modify LLM outputs or tool calls before continuing

{' ' * (self.list_depth - 1)}- [Interrupting tool calls](#interrupts-in-tools): Pause before executing tool calls to review and edit the tool call before execution

{' ' * (self.list_depth - 1)}- [Validating human input](#validating-human-input): Pause before proceeding to the next step to validate human input

### [​](#stream-with-human-in-the-loop-hitl-interrupts)Stream with human-in-the-loop (HITL) interrupts

When building interactive agents with human-in-the-loop workflows, you can stream both message chunks and node updates simultaneously to provide real-time feedback while handling interrupts.
Use multiple stream modes (`"messages"`, `"updates"`, and `"values"`) with `subgraphs=True` (if subgraphs are present) to:

{' ' * (self.list_depth - 1)}- Stream AI responses in real-time as they’re generated

{' ' * (self.list_depth - 1)}- Detect when the graph encounters an interrupt

{' ' * (self.list_depth - 1)}- Handle user input and resume execution seamlessly

```
from langchain.messages import AIMessageChunk
from langgraph.types import Command

for chunk in graph.stream(
 initial_input,
 stream_mode=["messages", "updates", "values"],
 subgraphs=True,
 config=config,
 version="v2",
):
 if chunk["type"] == "messages":
 msg, _ = chunk["data"]
 if isinstance(msg, AIMessageChunk) and msg.content:
 display_streaming_content(msg.content)

 elif chunk["type"] == "values" and chunk.get("interrupts"):
 interrupt_info = chunk["interrupts"][0].value
 user_response = get_user_input(interrupt_info)
 initial_input = Command(resume=user_response)
 break

 elif chunk["type"] == "updates":
 current_node = list(chunk["data"].keys())[0]

```

{' ' * (self.list_depth - 1)}- 
**`version="v2"`**: All chunks are `StreamPart` dicts with `type`, `ns`, and `data` keys; pending interrupts appear on `chunk["interrupts"]` for `"values"` parts

{' ' * (self.list_depth - 1)}- 
**`chunk["type"]`**: Narrow on the stream mode (`"messages"`, `"updates"`, `"values"`, etc.) for type inference

{' ' * (self.list_depth - 1)}- 
**`chunk["ns"]`**: Identifies the source graph (empty tuple for root, populated for subgraphs)

{' ' * (self.list_depth - 1)}- 
**`subgraphs=True`**: Required for interrupt detection in nested graphs

{' ' * (self.list_depth - 1)}- 
**`Command(resume=...)`**: Resumes graph execution with user-provided data

### [​](#handling-multiple-interrupts)Handling multiple interrupts

When parallel branches interrupt simultaneously (for example, fan-out to multiple nodes that each call `interrupt()`), you may need to resume multiple interrupts in a single invocation.
When resuming multiple interrupts with a single invocation, map each interrupt ID to its resume value.
This ensures each response is paired with the correct interrupt at runtime.

```
from typing import Annotated, TypedDict
import operator

from langgraph.checkpoint.memory import InMemorySaver
from langgraph.graph import END, START, StateGraph
from langgraph.types import Command, interrupt

class State(TypedDict):
 vals: Annotated[list[str], operator.add]

def node_a(state):
 answer = interrupt("question_a")
 return {"vals": [f"a:{answer}"]}

def node_b(state):
 answer = interrupt("question_b")
 return {"vals": [f"b:{answer}"]}

graph = (
 StateGraph(State)
 .add_node("a", node_a)
 .add_node("b", node_b)
 .add_edge(START, "a")
 .add_edge(START, "b")
 .add_edge("a", END)
 .add_edge("b", END)
 .compile(checkpointer=InMemorySaver())
)

config = {"configurable": {"thread_id": "1"}}

# Step 1: invoke with version="v2" to get typed GraphOutput
# Both parallel nodes hit interrupt() and pause
interrupted_result = graph.invoke({"vals": []}, config, version="v2")
print(interrupted_result.interrupts)
# > (Interrupt(value='question_a', id='...'), Interrupt(value='question_b', id='...'))

# Step 2: resume all pending interrupts at once
resume_map = {
 i.id: f"answer for {i.value}" for i in interrupted_result.interrupts
}
result = graph.invoke(Command(resume=resume_map), config, version="v2")

print("Final state:", result.value)
# Final state: {'vals': ['a:answer for question_a', 'b:answer for question_b']}

```

### [​](#approve-or-reject)Approve or reject

One of the most common uses of interrupts is to pause before a critical action and ask for approval. For example, you might want to ask a human to approve an API call, a database change, or any other important decision.

```
from typing import Literal
from langgraph.types import interrupt, Command

def approval_node(state: State) -> Command[Literal["proceed", "cancel"]]:
 # Pause execution; payload shows up in result.interrupts (v2) or result["__interrupt__"] (v1)
 is_approved = interrupt({
 "question": "Do you want to proceed with this action?",
 "details": state["action_details"]
 })

 # Route based on the response
 if is_approved:
 return Command(goto="proceed") # Runs after the resume payload is provided
 else:
 return Command(goto="cancel")

```

When you resume the graph, pass `True` to approve or `False` to reject:

```
# To approve
graph.invoke(Command(resume=True), config=config, version="v2")

# To reject
graph.invoke(Command(resume=False), config=config, version="v2")

```

Full example
```
from typing import Literal, Optional, TypedDict

from langgraph.checkpoint.memory import InMemorySaver
from langgraph.graph import END, START, StateGraph
from langgraph.types import Command, interrupt

class ApprovalState(TypedDict):
 action_details: str
 status: Optional[Literal["pending", "approved", "rejected"]]

def approval_node(state: ApprovalState) -> Command[Literal["proceed", "cancel"]]:
 # Expose details so the caller can render them in a UI
 decision = interrupt(
 {
 "question": "Approve this action?",
 "details": state["action_details"],
 }
 )

 # Route to the appropriate node after resume
 return Command(goto="proceed" if decision else "cancel")

def proceed_node(state: ApprovalState):
 return {"status": "approved"}

def cancel_node(state: ApprovalState):
 return {"status": "rejected"}

builder = StateGraph(ApprovalState)
builder.add_node("approval", approval_node)
builder.add_node("proceed", proceed_node)
builder.add_node("cancel", cancel_node)
builder.add_edge(START, "approval")
builder.add_edge("proceed", END)
builder.add_edge("cancel", END)

# Use a more durable checkpointer in production
checkpointer = InMemorySaver()
graph = builder.compile(checkpointer=checkpointer)

config = {"configurable": {"thread_id": "approval-123"}}
initial = graph.invoke(
 {"action_details": "Transfer $500", "status": "pending"},
 config=config,
 version="v2",
)
print(initial.interrupts) # -> (Interrupt(value={'question': ..., 'details': ...}),)

# Resume with the decision; True routes to proceed, False to cancel
resumed = graph.invoke(Command(resume=True), config=config, version="v2")
print(resumed.value["status"])

```

### [​](#review-and-edit-state)Review and edit state

Sometimes you want to let a human review and edit part of the graph state before continuing. This is useful for correcting LLMs, adding missing information, or making adjustments.

```
from langgraph.types import interrupt

def review_node(state: State):
 # Pause and show the current content for review (surfaces in result.interrupts with v2)
 edited_content = interrupt({
 "instruction": "Review and edit this content",
 "content": state["generated_text"]
 })

 # Update the state with the edited version
 return {"generated_text": edited_content}

```

When resuming, provide the edited content:

```
graph.invoke(
 Command(resume="The edited and improved text"), # Value becomes the return from interrupt()
 config=config,
 version="v2",
)

```

Full example
```
from typing import TypedDict

from langgraph.checkpoint.memory import MemorySaver
from langgraph.graph import END, START, StateGraph
from langgraph.types import Command, interrupt

class ReviewState(TypedDict):
 generated_text: str

def review_node(state: ReviewState):
 # Ask a reviewer to edit the generated content
 updated = interrupt(
 {
 "instruction": "Review and edit this content",
 "content": state["generated_text"],
 }
 )
 return {"generated_text": updated}

builder = StateGraph(ReviewState)
builder.add_node("review", review_node)
builder.add_edge(START, "review")
builder.add_edge("review", END)

checkpointer = MemorySaver()
graph = builder.compile(checkpointer=checkpointer)

config = {"configurable": {"thread_id": "review-42"}}
initial = graph.invoke({"generated_text": "Initial draft"}, config=config, version="v2")
print(initial.interrupts) # -> (Interrupt(value={'instruction': ..., 'content': ...}),)

# Resume with the edited text from the reviewer
final_state = graph.invoke(
 Command(resume="Improved draft after review"),
 config=config,
 version="v2",
)
print(final_state.value["generated_text"]) # -> "Improved draft after review"

```

### [​](#interrupts-in-tools)Interrupts in tools

You can also place interrupts directly inside tool functions. This makes the tool itself pause for approval whenever it’s called, and allows for human review and editing of the tool call before it is executed.
First, define a tool that uses [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt):

```
from langchain.tools import tool
from langgraph.types import interrupt

@tool
def send_email(to: str, subject: str, body: str):
 """Send an email to a recipient."""

 # Pause before sending; payload surfaces in result.interrupts (v2)
 response = interrupt({
 "action": "send_email",
 "to": to,
 "subject": subject,
 "body": body,
 "message": "Approve sending this email?"
 })

 if response.get("action") == "approve":
 # Resume value can override inputs before executing
 final_to = response.get("to", to)
 final_subject = response.get("subject", subject)
 final_body = response.get("body", body)
 return f"Email sent to {final_to} with subject '{final_subject}'"
 return "Email cancelled by user"

```

This approach is useful when you want the approval logic to live with the tool itself, making it reusable across different parts of your graph. The LLM can call the tool naturally, and the interrupt will pause execution whenever the tool is invoked, allowing you to approve, edit, or cancel the action.

Full example
```
import sqlite3
from typing import TypedDict

from langchain.tools import tool
from langchain_anthropic import ChatAnthropic
from langgraph.checkpoint.sqlite import SqliteSaver
from langgraph.graph import StateGraph, START, END
from langgraph.types import Command, interrupt

class AgentState(TypedDict):
 messages: list[dict]

@tool
def send_email(to: str, subject: str, body: str):
 """Send an email to a recipient."""

 # Pause before sending; payload surfaces in result.interrupts (v2)
 response = interrupt({
 "action": "send_email",
 "to": to,
 "subject": subject,
 "body": body,
 "message": "Approve sending this email?",
 })

 if response.get("action") == "approve":
 final_to = response.get("to", to)
 final_subject = response.get("subject", subject)
 final_body = response.get("body", body)

 # Actually send the email (your implementation here)
 print(f"[send_email] to={final_to} subject={final_subject} body={final_body}")
 return f"Email sent to {final_to}"

 return "Email cancelled by user"

model = ChatAnthropic(model="claude-sonnet-4-6").bind_tools([send_email])

def agent_node(state: AgentState):
 # LLM may decide to call the tool; interrupt pauses before sending
 result = model.invoke(state["messages"])
 return {"messages": state["messages"] + [result]}

builder = StateGraph(AgentState)
builder.add_node("agent", agent_node)
builder.add_edge(START, "agent")
builder.add_edge("agent", END)

checkpointer = SqliteSaver(sqlite3.connect("tool-approval.db"))
graph = builder.compile(checkpointer=checkpointer)

config = {"configurable": {"thread_id": "email-workflow"}}
initial = graph.invoke(
 {
 "messages": [
 {"role": "user", "content": "Send an email to alice@example.com about the meeting"}
 ]
 },
 config=config,
 version="v2",
)
print(initial.interrupts) # -> (Interrupt(value={'action': 'send_email', ...}),)

# Resume with approval and optionally edited arguments
resumed = graph.invoke(
 Command(resume={"action": "approve", "subject": "Updated subject"}),
 config=config,
 version="v2",
)
print(resumed.value["messages"][-1]) # -> Tool result returned by send_email

```

### [​](#validating-human-input)Validating human input

Sometimes you need to validate input from humans and ask again if it’s invalid. You can do this using multiple [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) calls in a loop.

```
from langgraph.types import interrupt

def get_age_node(state: State):
 prompt = "What is your age?"

 while True:
 answer = interrupt(prompt) # payload surfaces in result.interrupts (v2)

 # Validate the input
 if isinstance(answer, int) and answer > 0:
 # Valid input - continue
 break
 else:
 # Invalid input - ask again with a more specific prompt
 prompt = f"'{answer}' is not a valid age. Please enter a positive number."

 return {"age": answer}

```

Each time you resume the graph with invalid input, it will ask again with a clearer message. Once valid input is provided, the node completes and the graph continues.

Full example
```
import sqlite3
from typing import TypedDict

from langgraph.checkpoint.sqlite import SqliteSaver
from langgraph.graph import END, START, StateGraph
from langgraph.types import Command, interrupt

class FormState(TypedDict):
 age: int | None

def get_age_node(state: FormState):
 prompt = "What is your age?"

 while True:
 answer = interrupt(prompt)

 if isinstance(answer, int) and answer > 0:
 return {"age": answer}

 prompt = f"'{answer}' is not a valid age. Please enter a positive number."

builder = StateGraph(FormState)
builder.add_node("collect_age", get_age_node)
builder.add_edge(START, "collect_age")
builder.add_edge("collect_age", END)

checkpointer = SqliteSaver(sqlite3.connect("forms.db"))
graph = builder.compile(checkpointer=checkpointer)

config = {"configurable": {"thread_id": "form-1"}}
first = graph.invoke({"age": None}, config=config, version="v2")
print(first.interrupts) # -> (Interrupt(value='What is your age?', ...),)

# Provide invalid data; the node re-prompts
retry = graph.invoke(Command(resume="thirty"), config=config, version="v2")
print(retry.interrupts) # -> (Interrupt(value="'thirty' is not a valid age...", ...),)

# Provide valid data; loop exits and state updates
final = graph.invoke(Command(resume=30), config=config, version="v2")
print(final.value["age"]) # -> 30

```

## [​](#rules-of-interrupts)Rules of interrupts

When you call [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) within a node, LangGraph suspends execution by raising an exception that signals the runtime to pause. This exception propagates up through the call stack and is caught by the runtime, which notifies the graph to save the current state and wait for external input.
When execution resumes (after you provide the requested input), the runtime restarts the entire node from the beginning—it does not resume from the exact line where [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) was called. This means any code that ran before the [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) will execute again. Because of this, there’s a few important rules to follow when working with interrupts to ensure they behave as expected.

### [​](#do-not-wrap-interrupt-calls-in-try/except)Do not wrap `interrupt` calls in try/except

The way that [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) pauses execution at the point of the call is by throwing a special exception. If you wrap the [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) call in a try/except block, you will catch this exception and the interrupt will not be passed back to the graph.

{' ' * (self.list_depth - 1)}- ✅ Separate [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) calls from error-prone code

{' ' * (self.list_depth - 1)}- ✅ Use specific exception types in try/except blocks

Separating logicExplicit exception handling
```
def node_a(state: State):
 # ✅ Good: interrupting first, then handling
 # error conditions separately
 interrupt("What's your name?")
 try:
 fetch_data() # This can fail
 except Exception as e:
 print(e)
 return state

```

{' ' * (self.list_depth - 1)}- 🔴 Do not wrap [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) calls in bare try/except blocks

```
def node_a(state: State):
 # ❌ Bad: wrapping interrupt in bare try/except
 # will catch the interrupt exception
 try:
 interrupt("What's your name?")
 except Exception as e:
 print(e)
 return state

```

### [​](#do-not-reorder-interrupt-calls-within-a-node)Do not reorder `interrupt` calls within a node

It’s common to use multiple interrupts in a single node, however this can lead to unexpected behavior if not handled carefully.
When a node contains multiple interrupt calls, LangGraph keeps a list of resume values specific to the task executing the node. Whenever execution resumes, it starts at the beginning of the node. For each interrupt encountered, LangGraph checks if a matching value exists in the task’s resume list. Matching is **strictly index-based**, so the order of interrupt calls within the node is important.

{' ' * (self.list_depth - 1)}- ✅ Keep [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) calls consistent across node executions

```
def node_a(state: State):
 # ✅ Good: interrupt calls happen in the same order every time
 name = interrupt("What's your name?")
 age = interrupt("What's your age?")
 city = interrupt("What's your city?")

 return {
 "name": name,
 "age": age,
 "city": city
 }

```

{' ' * (self.list_depth - 1)}- 🔴 Do not conditionally skip [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) calls within a node

{' ' * (self.list_depth - 1)}- 🔴 Do not loop [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) calls using logic that isn’t deterministic across executions

Skipping interruptsLooping interrupts
```
def node_a(state: State):
 # ❌ Bad: conditionally skipping interrupts changes the order
 name = interrupt("What's your name?")

 # On first run, this might skip the interrupt
 # On resume, it might not skip it - causing index mismatch
 if state.get("needs_age"):
 age = interrupt("What's your age?")

 city = interrupt("What's your city?")

 return {"name": name, "city": city}

```

### [​](#do-not-return-complex-values-in-interrupt-calls)Do not return complex values in `interrupt` calls

Depending on which checkpointer is used, complex values may not be serializable (e.g. you can’t serialize a function). To make your graphs adaptable to any deployment, it’s best practice to only use values that can be reasonably serialized.

{' ' * (self.list_depth - 1)}- ✅ Pass simple, JSON-serializable types to [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt)

{' ' * (self.list_depth - 1)}- ✅ Pass dictionaries/objects with simple values

Simple valuesStructured data
```
def node_a(state: State):
 # ✅ Good: passing simple types that are serializable
 name = interrupt("What's your name?")
 count = interrupt(42)
 approved = interrupt(True)

 return {"name": name, "count": count, "approved": approved}

```

{' ' * (self.list_depth - 1)}- 🔴 Do not pass functions, class instances, or other complex objects to [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt)

FunctionsClass instances
```
def validate_input(value):
 return len(value) > 0

def node_a(state: State):
 # ❌ Bad: passing a function to interrupt
 # The function cannot be serialized
 response = interrupt({
 "question": "What's your name?",
 "validator": validate_input # This will fail
 })
 return {"name": response}

```

### [​](#side-effects-called-before-interrupt-must-be-idempotent)Side effects called before `interrupt` must be idempotent

Because interrupts work by re-running the nodes they were called from, side effects called before [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) should (ideally) be idempotent. For context, idempotency means that the same operation can be applied multiple times without changing the result beyond the initial execution.
As an example, you might have an API call to update a record inside of a node. If [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) is called after that call is made, it will be re-run multiple times when the node is resumed, potentially overwriting the initial update or creating duplicate records.

{' ' * (self.list_depth - 1)}- ✅ Use idempotent operations before [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt)

{' ' * (self.list_depth - 1)}- ✅ Place side effects after [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) calls

{' ' * (self.list_depth - 1)}- ✅ Separate side effects into separate nodes when possible

Idempotent operationsSide effects after interruptSeparating into different nodes
```
def node_a(state: State):
 # ✅ Good: using upsert operation which is idempotent
 # Running this multiple times will have the same result
 db.upsert_user(
 user_id=state["user_id"],
 status="pending_approval"
 )

 approved = interrupt("Approve this change?")

 return {"approved": approved}

```

{' ' * (self.list_depth - 1)}- 🔴 Do not perform non-idempotent operations before [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt)

{' ' * (self.list_depth - 1)}- 🔴 Do not create new records without checking if they exist

Creating recordsAppending to lists
```
def node_a(state: State):
 # ❌ Bad: creating a new record before interrupt
 # This will create duplicate records on each resume
 audit_id = db.create_audit_log({
 "user_id": state["user_id"],
 "action": "pending_approval",
 "timestamp": datetime.now()
 })

 approved = interrupt("Approve this change?")

 return {"approved": approved, "audit_id": audit_id}

```

## [​](#using-with-subgraphs-called-as-functions)Using with subgraphs called as functions

When invoking a subgraph within a node, the parent graph will resume execution from the **beginning of the node** where the subgraph was invoked and the [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) was triggered. Similarly, the **subgraph** will also resume from the beginning of the node where [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) was called.

```
def node_in_parent_graph(state: State):
 some_code() # <-- This will re-execute when resumed
 # Invoke a subgraph as a function.
 # The subgraph contains an `interrupt` call.
 subgraph_result = subgraph.invoke(some_input)
 # ...

def node_in_subgraph(state: State):
 some_other_code() # <-- This will also re-execute when resumed
 result = interrupt("What's your name?")
 # ...

```

## [​](#debugging-with-interrupts)Debugging with interrupts

To debug and test a graph, you can use static interrupts as breakpoints to step through the graph execution one node at a time. Static interrupts are triggered at defined points either before or after a node executes. You can set these by specifying `interrupt_before` and `interrupt_after` when compiling the graph.
Static interrupts are **not** recommended for human-in-the-loop workflows. Use the [`interrupt`](https://reference.langchain.com/python/langgraph/types/interrupt) function instead.

{' ' * (self.list_depth - 1)}- At compile time
{' ' * (self.list_depth - 1)}- At run time
```
graph = builder.compile(
 interrupt_before=["node_a"],
 interrupt_after=["node_b", "node_c"],
 checkpointer=checkpointer,
)

# Pass a thread ID to the graph
config = {
 "configurable": {
 "thread_id": "some_thread"
 }
}

# Run the graph until the breakpoint
graph.invoke(inputs, config=config)

# Resume the graph
graph.invoke(None, config=config)

```

{' ' * (self.list_depth - 1)}- The breakpoints are set during `compile` time.

{' ' * (self.list_depth - 1)}- `interrupt_before` specifies the nodes where execution should pause before the node is executed.

{' ' * (self.list_depth - 1)}- `interrupt_after` specifies the nodes where execution should pause after the node is executed.

{' ' * (self.list_depth - 1)}- A checkpointer is required to enable breakpoints.

{' ' * (self.list_depth - 1)}- The graph is run until the first breakpoint is hit.

{' ' * (self.list_depth - 1)}- The graph is resumed by passing in `None` for the input. This will run the graph until the next breakpoint is hit.

```
config = {
 "configurable": {
 "thread_id": "some_thread"
 }
}

# Run the graph until the breakpoint
graph.invoke(
 inputs,
 interrupt_before=["node_a"],
 interrupt_after=["node_b", "node_c"],
 config=config,
)

# Resume the graph
graph.invoke(None, config=config)

```

{' ' * (self.list_depth - 1)}- `graph.invoke` is called with the `interrupt_before` and `interrupt_after` parameters. This is a run-time configuration and can be changed for every invocation.

{' ' * (self.list_depth - 1)}- `interrupt_before` specifies the nodes where execution should pause before the node is executed.

{' ' * (self.list_depth - 1)}- `interrupt_after` specifies the nodes where execution should pause after the node is executed.

{' ' * (self.list_depth - 1)}- The graph is run until the first breakpoint is hit.

{' ' * (self.list_depth - 1)}- The graph is resumed by passing in `None` for the input. This will run the graph until the next breakpoint is hit.

To debug your interrupts, use [LangSmith](/langsmith/home).

### [​](#using-langsmith-studio)Using LangSmith Studio

You can use [LangSmith Studio](/langsmith/studio) to set static interrupts in your graph in the UI before running the graph. You can also use the UI to inspect the graph state at any point in the execution.

![image](https://mintcdn.com/langchain-5e9cc07a/dL5Sn6Cmy9pwtY0V/oss/images/static-interrupt.png?fit=max&auto=format&n=dL5Sn6Cmy9pwtY0V&q=85&s=5aa4e7cea2ab147cef5b4e210dd6c4a1)

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langgraph/interrupts.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[StreamingPrevious](/oss/python/langgraph/streaming)[Use time-travelNext](/oss/python/langgraph/use-time-travel)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
