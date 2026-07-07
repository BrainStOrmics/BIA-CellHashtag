# Build a personal assistant with subagents

[Tutorials](/oss/python/deepagents/data-analysis)[Multi-agent](/oss/python/langchain/multi-agent/subagents-personal-assistant)

# Build a personal assistant with subagents
Copy pageCopy page
> 

## Documentation Index

Fetch the complete documentation index at: [https://docs.langchain.com/llms.txt](https://docs.langchain.com/llms.txt)

Use this file to discover all available pages before exploring further.

## [​](#overview)Overview

The **supervisor pattern** is a [multi-agent](/oss/python/langchain/multi-agent) architecture where a central supervisor agent coordinates specialized worker agents. This approach excels when tasks require different types of expertise. Rather than building one agent that manages tool selection across domains, you create focused specialists coordinated by a supervisor who understands the overall workflow.
In this tutorial, you’ll build a personal assistant system that demonstrates these benefits through a realistic workflow. The system will coordinate two specialists with fundamentally different responsibilities:

{' ' * (self.list_depth - 1)}- A **calendar agent** that handles scheduling, availability checking, and event management.

{' ' * (self.list_depth - 1)}- An **email agent** that manages communication, drafts messages, and sends notifications.

We will also incorporate [human-in-the-loop review](/oss/python/langchain/human-in-the-loop) to allow users to approve, edit, and reject actions (such as outbound emails) as desired.

### [​](#why-use-a-supervisor)Why use a supervisor?

Multi-agent architectures allow you to partition [tools](/oss/python/langchain/tools) across workers, each with their own individual prompts or instructions. Consider an agent with direct access to all calendar and email APIs: it must choose from many similar tools, understand exact formats for each API, and handle multiple domains simultaneously. If performance degrades, it may be helpful to separate related tools and associated prompts into logical groups (in part to manage iterative improvements).

### [​](#concepts)Concepts

We will cover the following concepts:

{' ' * (self.list_depth - 1)}- [Multi-agent systems](/oss/python/langchain/multi-agent)

{' ' * (self.list_depth - 1)}- [Human-in-the-loop review](/oss/python/langchain/human-in-the-loop)

## [​](#setup)Setup

### [​](#installation)Installation

This tutorial requires the `langchain` package:
pipconda
```
pip install langchain

```

For more details, see our [Installation guide](/oss/python/langchain/install).

### [​](#langsmith)LangSmith

Set up [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-multi-agent-subagents-personal-assistant) to inspect what is happening inside your agent. Then set the following environment variables:
bashpython
```
export LANGSMITH_TRACING="true"
export LANGSMITH_API_KEY="..."

```

### [​](#components)Components

We will need to select a chat model from LangChain’s suite of integrations:

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

## [​](#1-define-tools)1. Define tools

Start by defining the tools that require structured inputs. In real applications, these would call actual APIs (Google Calendar, SendGrid, etc.). For this tutorial, you’ll use stubs to demonstrate the pattern.

```
from langchain.tools import tool

@tool
def create_calendar_event(
 title: str,
 start_time: str, # ISO format: "2024-01-15T14:00:00"
 end_time: str, # ISO format: "2024-01-15T15:00:00"
 attendees: list[str], # email addresses
 location: str = ""
) -> str:
 """Create a calendar event. Requires exact ISO datetime format."""
 # Stub: In practice, this would call Google Calendar API, Outlook API, etc.
 return f"Event created: {title} from {start_time} to {end_time} with {len(attendees)} attendees"

@tool
def send_email(
 to: list[str], # email addresses
 subject: str,
 body: str,
 cc: list[str] = []
) -> str:
 """Send an email via email API. Requires properly formatted addresses."""
 # Stub: In practice, this would call SendGrid, Gmail API, etc.
 return f"Email sent to {', '.join(to)} - Subject: {subject}"

@tool
def get_available_time_slots(
 attendees: list[str],
 date: str, # ISO format: "2024-01-15"
 duration_minutes: int
) -> list[str]:
 """Check calendar availability for given attendees on a specific date."""
 # Stub: In practice, this would query calendar APIs
 return ["09:00", "14:00", "16:00"]

```

## [​](#2-create-specialized-sub-agents)2. Create specialized sub-agents

Next, we’ll create specialized sub-agents that handle each domain.

### [​](#create-a-calendar-agent)Create a calendar agent

The calendar agent understands natural language scheduling requests and translates them into precise API calls. It handles date parsing, availability checking, and event creation.

```
from langchain.agents import create_agent

CALENDAR_AGENT_PROMPT = (
 "You are a calendar scheduling assistant. "
 "Parse natural language scheduling requests (e.g., 'next Tuesday at 2pm') "
 "into proper ISO datetime formats. "
 "Use get_available_time_slots to check availability when needed. "
 "If there is no suitable time slot, stop and confirm unavailability in your response. "
 "Use create_calendar_event to schedule events. "
 "Always confirm what was scheduled in your final response."
)

calendar_agent = create_agent(
 model,
 tools=[create_calendar_event, get_available_time_slots],
 system_prompt=CALENDAR_AGENT_PROMPT,
)

```

Test the calendar agent to see how it handles natural language scheduling:

```
query = "Schedule a team meeting next Tuesday at 2pm for 1 hour"

for step in calendar_agent.stream(
 {"messages": [{"role": "user", "content": query}]}
):
 for update in step.values():
 for message in update.get("messages", []):
 message.pretty_print()

```

```
================================== Ai Message ==================================
Tool Calls:
 get_available_time_slots (call_EIeoeIi1hE2VmwZSfHStGmXp)
 Call ID: call_EIeoeIi1hE2VmwZSfHStGmXp
 Args:
 attendees: []
 date: 2024-06-18
 duration_minutes: 60
================================= Tool Message =================================
Name: get_available_time_slots

["09:00", "14:00", "16:00"]
================================== Ai Message ==================================
Tool Calls:
 create_calendar_event (call_zgx3iJA66Ut0W8S3NpT93kEB)
 Call ID: call_zgx3iJA66Ut0W8S3NpT93kEB
 Args:
 title: Team Meeting
 start_time: 2024-06-18T14:00:00
 end_time: 2024-06-18T15:00:00
 attendees: []
================================= Tool Message =================================
Name: create_calendar_event

Event created: Team Meeting from 2024-06-18T14:00:00 to 2024-06-18T15:00:00 with 0 attendees
================================== Ai Message ==================================

The team meeting has been scheduled for next Tuesday, June 18th, at 2:00 PM and will last for 1 hour. If you need to add attendees or a location, please let me know!

```

The agent parses “next Tuesday at 2pm” into ISO format (“2024-01-16T14:00:00”), calculates the end time, calls `create_calendar_event`, and returns a natural language confirmation.

### [​](#create-an-email-agent)Create an email agent

The email agent handles message composition and sending. It focuses on extracting recipient information, crafting appropriate subject lines and body text, and managing email communication.

```
EMAIL_AGENT_PROMPT = (
 "You are an email assistant. "
 "Compose professional emails based on natural language requests. "
 "Extract recipient information and craft appropriate subject lines and body text. "
 "Use send_email to send the message. "
 "Always confirm what was sent in your final response."
)

email_agent = create_agent(
 model,
 tools=[send_email],
 system_prompt=EMAIL_AGENT_PROMPT,
)

```

Test the email agent with a natural language request:

```
query = "Send the design team a reminder about reviewing the new mockups"

for step in email_agent.stream(
 {"messages": [{"role": "user", "content": query}]}
):
 for update in step.values():
 for message in update.get("messages", []):
 message.pretty_print()

```

```
================================== Ai Message ==================================
Tool Calls:
 send_email (call_OMl51FziTVY6CRZvzYfjYOZr)
 Call ID: call_OMl51FziTVY6CRZvzYfjYOZr
 Args:
 to: ['design-team@example.com']
 subject: Reminder: Please Review the New Mockups
 body: Hi Design Team,

This is a friendly reminder to review the new mockups at your earliest convenience. Your feedback is important to ensure that we stay on track with our project timeline.

Please let me know if you have any questions or need additional information.

Thank you!

Best regards,
================================= Tool Message =================================
Name: send_email

Email sent to design-team@example.com - Subject: Reminder: Please Review the New Mockups
================================== Ai Message ==================================

I've sent a reminder to the design team asking them to review the new mockups. If you need any further communication on this topic, just let me know!

```

The agent infers the recipient from the informal request, crafts a professional subject line and body, calls `send_email`, and returns a confirmation. Each sub-agent has a narrow focus with domain-specific tools and prompts, allowing it to excel at its specific task.

## [​](#3-wrap-sub-agents-as-tools)3. Wrap sub-agents as tools

Now wrap each sub-agent as a tool that the supervisor can invoke. This is the key architectural step that creates the layered system. The supervisor will see high-level tools like “schedule_event”, not low-level tools like “create_calendar_event”.

```
@tool
def schedule_event(request: str) -> str:
 """Schedule calendar events using natural language.

 Use this when the user wants to create, modify, or check calendar appointments.
 Handles date/time parsing, availability checking, and event creation.

 Input: Natural language scheduling request (e.g., 'meeting with design team
 next Tuesday at 2pm')
 """
 result = calendar_agent.invoke({
 "messages": [{"role": "user", "content": request}]
 })
 return result["messages"][-1].text

@tool
def manage_email(request: str) -> str:
 """Send emails using natural language.

 Use this when the user wants to send notifications, reminders, or any email
 communication. Handles recipient extraction, subject generation, and email
 composition.

 Input: Natural language email request (e.g., 'send them a reminder about
 the meeting')
 """
 result = email_agent.invoke({
 "messages": [{"role": "user", "content": request}]
 })
 return result["messages"][-1].text

```

The tool descriptions help the supervisor decide when to use each tool, so make them clear and specific. We return only the sub-agent’s final response, as the supervisor doesn’t need to see intermediate reasoning or tool calls.

## [​](#4-create-the-supervisor-agent)4. Create the supervisor agent

Now create the supervisor that orchestrates the sub-agents. The supervisor only sees high-level tools and makes routing decisions at the domain level, not the individual API level.

```
SUPERVISOR_PROMPT = (
 "You are a helpful personal assistant. "
 "You can schedule calendar events and send emails. "
 "Break down user requests into appropriate tool calls and coordinate the results. "
 "When a request involves multiple actions, use multiple tools in sequence."
)

supervisor_agent = create_agent(
 model,
 tools=[schedule_event, manage_email],
 system_prompt=SUPERVISOR_PROMPT,
)

```

## [​](#5-use-the-supervisor)5. Use the supervisor

Now test your complete system with complex requests that require coordination across multiple domains:

### [​](#example-1-simple-single-domain-request)Example 1: Simple single-domain request

```
query = "Schedule a team standup for tomorrow at 9am"

for step in supervisor_agent.stream(
 {"messages": [{"role": "user", "content": query}]}
):
 for update in step.values():
 for message in update.get("messages", []):
 message.pretty_print()

```

```
================================== Ai Message ==================================
Tool Calls:
 schedule_event (call_mXFJJDU8bKZadNUZPaag8Lct)
 Call ID: call_mXFJJDU8bKZadNUZPaag8Lct
 Args:
 request: Schedule a team standup for tomorrow at 9am with Alice and Bob.
================================= Tool Message =================================
Name: schedule_event

The team standup has been scheduled for tomorrow at 9:00 AM with Alice and Bob. If you need to make any changes or add more details, just let me know!
================================== Ai Message ==================================

The team standup with Alice and Bob is scheduled for tomorrow at 9:00 AM. If you need any further arrangements or adjustments, please let me know!

```

The supervisor identifies this as a calendar task, calls `schedule_event`, and the calendar agent handles date parsing and event creation.
For full transparency into the information flow, including prompts and responses for each chat model call, check out the [LangSmith trace](https://smith.langchain.com/public/91a9a95f-fba9-4e84-aff0-371861ad2f4a/r) for the above run.

### [​](#example-2-complex-multi-domain-request)Example 2: Complex multi-domain request

```
query = (
 "Schedule a meeting with the design team next Tuesday at 2pm for 1 hour, "
 "and send them an email reminder about reviewing the new mockups."
)

for step in supervisor_agent.stream(
 {"messages": [{"role": "user", "content": query}]}
):
 for update in step.values():
 for message in update.get("messages", []):
 message.pretty_print()

```

```
================================== Ai Message ==================================
Tool Calls:
 schedule_event (call_YA68mqF0koZItCFPx0kGQfZi)
 Call ID: call_YA68mqF0koZItCFPx0kGQfZi
 Args:
 request: meeting with the design team next Tuesday at 2pm for 1 hour
 manage_email (call_XxqcJBvVIuKuRK794ZIzlLxx)
 Call ID: call_XxqcJBvVIuKuRK794ZIzlLxx
 Args:
 request: send the design team an email reminder about reviewing the new mockups
================================= Tool Message =================================
Name: schedule_event

Your meeting with the design team is scheduled for next Tuesday, June 18th, from 2:00pm to 3:00pm. Let me know if you need to add more details or make any changes!
================================= Tool Message =================================
Name: manage_email

I've sent an email reminder to the design team requesting them to review the new mockups. If you need to include more information or recipients, just let me know!
================================== Ai Message ==================================

Your meeting with the design team is scheduled for next Tuesday, June 18th, from 2:00pm to 3:00pm.

I've also sent an email reminder to the design team, asking them to review the new mockups.

Let me know if you'd like to add more details to the meeting or include additional information in the email!

```

The supervisor recognizes this requires both calendar and email actions, calls `schedule_event` for the meeting, then calls `manage_email` for the reminder. Each sub-agent completes its task, and the supervisor synthesizes both results into a coherent response.
Refer to the [LangSmith trace](https://smith.langchain.com/public/95cd00a3-d1f9-4dba-9731-7bf733fb6a3c/r) to see the detailed information flow for the above run, including individual chat model prompts and responses.

### [​](#complete-working-example)Complete working example

Here’s everything together in a runnable script:

Show View complete code
```
"""
Personal Assistant Supervisor Example

This example demonstrates the tool calling pattern for multi-agent systems.
A supervisor agent coordinates specialized sub-agents (calendar and email)
that are wrapped as tools.
"""

from langchain.tools import tool
from langchain.agents import create_agent
from langchain.chat_models import init_chat_model

# ============================================================================
# Step 1: Define low-level API tools (stubbed)
# ============================================================================

@tool
def create_calendar_event(
 title: str,
 start_time: str, # ISO format: "2024-01-15T14:00:00"
 end_time: str, # ISO format: "2024-01-15T15:00:00"
 attendees: list[str], # email addresses
 location: str = ""
) -> str:
 """Create a calendar event. Requires exact ISO datetime format."""
 return f"Event created: {title} from {start_time} to {end_time} with {len(attendees)} attendees"

@tool
def send_email(
 to: list[str], # email addresses
 subject: str,
 body: str,
 cc: list[str] = []
) -> str:
 """Send an email via email API. Requires properly formatted addresses."""
 return f"Email sent to {', '.join(to)} - Subject: {subject}"

@tool
def get_available_time_slots(
 attendees: list[str],
 date: str, # ISO format: "2024-01-15"
 duration_minutes: int
) -> list[str]:
 """Check calendar availability for given attendees on a specific date."""
 return ["09:00", "14:00", "16:00"]

# ============================================================================
# Step 2: Create specialized sub-agents
# ============================================================================

model = init_chat_model("gpt-5.4") # for example

calendar_agent = create_agent(
 model,
 tools=[create_calendar_event, get_available_time_slots],
 system_prompt=(
 "You are a calendar scheduling assistant. "
 "Parse natural language scheduling requests (e.g., 'next Tuesday at 2pm') "
 "into proper ISO datetime formats. "
 "Use get_available_time_slots to check availability when needed. "
 "If there is no suitable time slot, stop and confirm unavailability in your response. "
 "Use create_calendar_event to schedule events. "
 "Always confirm what was scheduled in your final response."
 )
)

email_agent = create_agent(
 model,
 tools=[send_email],
 system_prompt=(
 "You are an email assistant. "
 "Compose professional emails based on natural language requests. "
 "Extract recipient information and craft appropriate subject lines and body text. "
 "Use send_email to send the message. "
 "Always confirm what was sent in your final response."
 )
)

# ============================================================================
# Step 3: Wrap sub-agents as tools for the supervisor
# ============================================================================

@tool
def schedule_event(request: str) -> str:
 """Schedule calendar events using natural language.

 Use this when the user wants to create, modify, or check calendar appointments.
 Handles date/time parsing, availability checking, and event creation.

 Input: Natural language scheduling request (e.g., 'meeting with design team
 next Tuesday at 2pm')
 """
 result = calendar_agent.invoke({
 "messages": [{"role": "user", "content": request}]
 })
 return result["messages"][-1].text

@tool
def manage_email(request: str) -> str:
 """Send emails using natural language.

 Use this when the user wants to send notifications, reminders, or any email
 communication. Handles recipient extraction, subject generation, and email
 composition.

 Input: Natural language email request (e.g., 'send them a reminder about
 the meeting')
 """
 result = email_agent.invoke({
 "messages": [{"role": "user", "content": request}]
 })
 return result["messages"][-1].text

# ============================================================================
# Step 4: Create the supervisor agent
# ============================================================================

supervisor_agent = create_agent(
 model,
 tools=[schedule_event, manage_email],
 system_prompt=(
 "You are a helpful personal assistant. "
 "You can schedule calendar events and send emails. "
 "Break down user requests into appropriate tool calls and coordinate the results. "
 "When a request involves multiple actions, use multiple tools in sequence."
 )
)

# ============================================================================
# Step 5: Use the supervisor
# ============================================================================

if __name__ == "__main__":
 # Example: User request requiring both calendar and email coordination
 user_request = (
 "Schedule a meeting with the design team next Tuesday at 2pm for 1 hour, "
 "and send them an email reminder about reviewing the new mockups."
 )

 print("User Request:", user_request)
 print("\n" + "="*80 + "\n")

 for step in supervisor_agent.stream(
 {"messages": [{"role": "user", "content": user_request}]}
 ):
 for update in step.values():
 for message in update.get("messages", []):
 message.pretty_print()

```

### [​](#understanding-the-architecture)Understanding the architecture

Your system has three layers. The bottom layer contains rigid API tools that require exact formats. The middle layer contains sub-agents that accept natural language, translate it to structured API calls, and return natural language confirmations. The top layer contains the supervisor that routes to high-level capabilities and synthesizes results.
This separation of concerns provides several benefits: each layer has a focused responsibility, you can add new domains without affecting existing ones, and you can test and iterate on each layer independently.

## [​](#6-add-human-in-the-loop-review)6. Add human-in-the-loop review

It can be prudent to incorporate [human-in-the-loop review](/oss/python/langchain/human-in-the-loop) of sensitive actions. LangChain includes [built-in middleware](/oss/python/langchain/human-in-the-loop#configuring-interrupts) to review tool calls, in this case the tools invoked by sub-agents.
Let’s add human-in-the-loop review to both sub-agents:

{' ' * (self.list_depth - 1)}- We configure the `create_calendar_event` and `send_email` tools to interrupt, permitting all [response types](/oss/python/langchain/human-in-the-loop) (`approve`, `edit`, `reject`)

{' ' * (self.list_depth - 1)}- We add a [checkpointer](/oss/python/langchain/short-term-memory) **only to the top-level agent**. This is required to pause and resume execution.

```
from langchain.agents import create_agent
from langchain.agents.middleware import HumanInTheLoopMiddleware 
from langgraph.checkpoint.memory import InMemorySaver 

calendar_agent = create_agent(
 model,
 tools=[create_calendar_event, get_available_time_slots],
 system_prompt=CALENDAR_AGENT_PROMPT,
 middleware=[
 HumanInTheLoopMiddleware(
 interrupt_on={"create_calendar_event": True},
 description_prefix="Calendar event pending approval",
 ),
 ],
)

email_agent = create_agent(
 model,
 tools=[send_email],
 system_prompt=EMAIL_AGENT_PROMPT,
 middleware=[
 HumanInTheLoopMiddleware(
 interrupt_on={"send_email": True},
 description_prefix="Outbound email pending approval",
 ),
 ],
)

supervisor_agent = create_agent(
 model,
 tools=[schedule_event, manage_email],
 system_prompt=SUPERVISOR_PROMPT,
 checkpointer=InMemorySaver(),
)

```

Let’s repeat the query. Note that we gather interrupt events into a list to access downstream:

```
query = (
 "Schedule a meeting with the design team next Tuesday at 2pm for 1 hour, "
 "and send them an email reminder about reviewing the new mockups."
)

config = {"configurable": {"thread_id": "6"}}

interrupts = []
for step in supervisor_agent.stream(
 {"messages": [{"role": "user", "content": query}]},
 config,
):
 for update in step.values():
 if isinstance(update, dict):
 for message in update.get("messages", []):
 message.pretty_print()
 else:
 interrupt_ = update[0]
 interrupts.append(interrupt_)
 print(f"\nINTERRUPTED: {interrupt_.id}")

```

```
================================== Ai Message ==================================
Tool Calls:
 schedule_event (call_t4Wyn32ohaShpEZKuzZbl83z)
 Call ID: call_t4Wyn32ohaShpEZKuzZbl83z
 Args:
 request: Schedule a meeting with the design team next Tuesday at 2pm for 1 hour.
 manage_email (call_JWj4vDJ5VMnvkySymhCBm4IR)
 Call ID: call_JWj4vDJ5VMnvkySymhCBm4IR
 Args:
 request: Send an email reminder to the design team about reviewing the new mockups before our meeting next Tuesday at 2pm.

INTERRUPTED: 4f994c9721682a292af303ec1a46abb7

INTERRUPTED: 2b56f299be313ad8bc689eff02973f16

```

This time we’ve interrupted execution. Let’s inspect the interrupt events:

```
for interrupt_ in interrupts:
 for request in interrupt_.value["action_requests"]:
 print(f"INTERRUPTED: {interrupt_.id}")
 print(f"{request['description']}\n")

```

```
INTERRUPTED: 4f994c9721682a292af303ec1a46abb7
Calendar event pending approval

Tool: create_calendar_event
Args: {'title': 'Meeting with the Design Team', 'start_time': '2024-06-18T14:00:00', 'end_time': '2024-06-18T15:00:00', 'attendees': ['design team']}

INTERRUPTED: 2b56f299be313ad8bc689eff02973f16
Outbound email pending approval

Tool: send_email
Args: {'to': ['designteam@example.com'], 'subject': 'Reminder: Review New Mockups Before Meeting Next Tuesday at 2pm', 'body': "Hello Team,\n\nThis is a reminder to review the new mockups ahead of our meeting scheduled for next Tuesday at 2pm. Your feedback and insights will be valuable for our discussion and next steps.\n\nPlease ensure you've gone through the designs and are ready to share your thoughts during the meeting.\n\nThank you!\n\nBest regards,\n[Your Name]"}

```

We can specify decisions for each interrupt by referring to its ID using a [`Command`](https://reference.langchain.com/python/langgraph/types/Command). Refer to the [human-in-the-loop guide](/oss/python/langchain/human-in-the-loop) for additional details. For demonstration purposes, here we will accept the calendar event, but edit the subject of the outbound email:

```
from langgraph.types import Command 

resume = {}
for interrupt_ in interrupts:
 if interrupt_.id == "2b56f299be313ad8bc689eff02973f16":
 # Edit email
 edited_action = interrupt_.value["action_requests"][0].copy()
 edited_action["args"]["subject"] = "Mockups reminder"
 resume[interrupt_.id] = {
 "decisions": [{"type": "edit", "edited_action": edited_action}]
 }
 else:
 resume[interrupt_.id] = {"decisions": [{"type": "approve"}]}

interrupts = []
for step in supervisor_agent.stream(
 Command(resume=resume),
 config,
):
 for update in step.values():
 if isinstance(update, dict):
 for message in update.get("messages", []):
 message.pretty_print()
 else:
 interrupt_ = update[0]
 interrupts.append(interrupt_)
 print(f"\nINTERRUPTED: {interrupt_.id}")

```

```
================================= Tool Message =================================
Name: schedule_event

Your meeting with the design team has been scheduled for next Tuesday, June 18th, from 2:00 pm to 3:00 pm.
================================= Tool Message =================================
Name: manage_email

Your email reminder to the design team has been sent. Here’s what was sent:

- Recipient: designteam@example.com
- Subject: Mockups reminder
- Body: A reminder to review the new mockups before the meeting next Tuesday at 2pm, with a request for feedback and readiness for discussion.

Let me know if you need any further assistance!
================================== Ai Message ==================================

- Your meeting with the design team has been scheduled for next Tuesday, June 18th, from 2:00 pm to 3:00 pm.
- An email reminder has been sent to the design team about reviewing the new mockups before the meeting.

Let me know if you need any further assistance!

```

The run proceeds with our input.

## [​](#7-advanced-control-information-flow)7. Advanced: Control information flow

By default, sub-agents receive only the request string from the supervisor. You might want to pass additional context, such as conversation history or user preferences.

### [​](#pass-additional-conversational-context-to-sub-agents)Pass additional conversational context to sub-agents

```
from langchain.tools import tool, ToolRuntime

@tool
def schedule_event(
 request: str,
 runtime: ToolRuntime
) -> str:
 """Schedule calendar events using natural language."""
 # Customize context received by sub-agent
 original_user_message = next(
 message for message in runtime.state["messages"]
 if message.type == "human"
 )
 prompt = (
 "You are assisting with the following user inquiry:\n\n"
 f"{original_user_message.text}\n\n"
 "You are tasked with the following sub-request:\n\n"
 f"{request}"
 )
 result = calendar_agent.invoke({
 "messages": [{"role": "user", "content": prompt}],
 })
 return result["messages"][-1].text

```

This allows sub-agents to see the full conversation context, which can be useful for resolving ambiguities like “schedule it for the same time tomorrow” (referencing a previous conversation).
You can see the full context received by the sub agent in the [chat model call](https://smith.langchain.com/public/c7d54882-afb8-4039-9c5a-4112d0f458b0/r/6803571e-af78-4c68-904a-ecf55771084d) of the LangSmith trace.

### [​](#control-what-supervisor-receives)Control what supervisor receives

You can also customize what information flows back to the supervisor:

```
import json

@tool
def schedule_event(request: str) -> str:
 """Schedule calendar events using natural language."""
 result = calendar_agent.invoke({
 "messages": [{"role": "user", "content": request}]
 })

 # Option 1: Return just the confirmation message
 return result["messages"][-1].text

 # Option 2: Return structured data
 # return json.dumps({
 # "status": "success",
 # "event_id": "evt_123",
 # "summary": result["messages"][-1].text
 # })

```

**Important:** Make sure sub-agent prompts emphasize that their final message should contain all relevant information. A common failure mode is sub-agents that perform tool calls but don’t include the results in their final response.

## [​](#8-key-takeaways)8. Key takeaways

The supervisor pattern creates layers of abstraction where each layer has a clear responsibility. When designing a supervisor system, start with clear domain boundaries and give each sub-agent focused tools and prompts. Write clear tool descriptions for the supervisor, test each layer independently before integration, and control information flow based on your specific needs.
**When to use the supervisor pattern**Use the supervisor pattern when you have multiple distinct domains (calendar, email, CRM, database), each domain has multiple tools or complex logic, you want centralized workflow control, and sub-agents don’t need to converse directly with users.For simpler cases with just a few tools, use a single agent. When agents need to have conversations with users, use [handoffs](/oss/python/langchain/multi-agent/handoffs) instead. For peer-to-peer collaboration between agents, consider other multi-agent patterns.

## [​](#next-steps)Next steps

Learn about [handoffs](/oss/python/langchain/multi-agent/handoffs) for agent-to-agent conversations, explore [context engineering](/oss/python/langchain/context-engineering) to fine-tune information flow, read the [multi-agent overview](/oss/python/langchain/multi-agent) to compare different patterns, and use [LangSmith](https://smith.langchain.com?utm_source=docs&utm_medium=cta&utm_campaign=langsmith-signup&utm_content=oss-langchain-multi-agent-subagents-personal-assistant) to debug and monitor your multi-agent system.

---

[Connect these docs](/use-these-docs) to Claude, VSCode, and more via MCP for real-time answers.[Edit this page on GitHub](https://github.com/langchain-ai/docs/edit/main/src/oss/langchain/multi-agent/subagents-personal-assistant.mdx) or [file an issue](https://github.com/langchain-ai/docs/issues/new/choose).

Was this page helpful?YesNo[Build a voice agent with LangChainPrevious](/oss/python/langchain/voice-agent)[Build customer support with handoffsNext](/oss/python/langchain/multi-agent/handoffs-customer-support)⌘I[Docs by LangChain home page
![light logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-dark-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=5babf1a1962208fd7eed942fa2432ecb)
![dark logo](https://mintcdn.com/langchain-5e9cc07a/nQm-sjd_MByLhgeW/images/brand/langchain-docs-light-blue.png?fit=max&auto=format&n=nQm-sjd_MByLhgeW&q=85&s=0bcd2a1f2599ed228bcedf0f535b45b1)](/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)

Resources[Forum](https://forum.langchain.com/)[Changelog](https://changelog.langchain.com/)[LangChain Academy](https://academy.langchain.com/)[Contact Sales](https://www.langchain.com/contact-sales)

Company[Home](https://langchain.com/)[Trust Center](https://trust.langchain.com/)[Careers](https://langchain.com/careers)[Blog](https://blog.langchain.com/)[github](https://github.com/langchain-ai)[x](https://x.com/LangChain)[linkedin](https://www.linkedin.com/company/langchain)[youtube](https://www.youtube.com/@LangChain)
